import json
import time
import traceback
from pathlib import Path

from tqdm.auto import tqdm
from metadata import Metadata
from reader import Reader
from writer import Writer
from api import API
from model import XGBoostModel
from database import Database
import argparse
from datetime import datetime


METADATA_FILENAME = "ena_metadata.tsv"

RANKS = [
    "biotype",
    "clade",
    "class",
    "family",
    "genus",
    "kingdom",
    "no rank",
    "order",
    "pathogroup",
    "phylum",
    "serogroup",
    "serotype",
    "species",
    "species group",
    "species subgroup",
    "strain",
    "subclass",
    "subfamily",
    "subgenus",
    "suborder",
    "subspecies",
    "superkingdom",
    "tribe",
]


def format_duration(seconds):
    """Convert seconds to a formatted duration string (hh:mm:ss)."""
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}"


def create_db(conf: dict):
    input_path = (Path(conf["allthebacteria"]["assembly_dir"]),)
    output_path = (Path(conf["db"]["output_dir"]),)
    metadata_path = (Path(conf["allthebacteria"]["metadata_dir"]) / METADATA_FILENAME,)
    api_cache_path = Path(conf["api"]["cache_dir"])
    output_path.mkdir(parents=True, exist_ok=True)

    json_log_path = output_path / "error_log.json"

    if json_log_path.exists():
        with open(json_log_path, "r", encoding="utf-8") as error_log_file:
            error_log = json.load(error_log_file)
    else:
        error_log = {"complete": [], "incomplete": [], "duration": {}}

    archives = list(input_path.glob("*.xz"))

    api = API(
        api_cache_dir=api_cache_path,
        email=conf["api"]["email"],
        can_download=conf["api"]["can_download"],
    )
    md = Metadata(csv_path=metadata_path, api=api)
    md.md  # preload
    reader = Reader(md)
    writer = Writer(output_path)

    for archive_path in tqdm(
        archives, desc=f"Processing {str(input_path)} -> {str(output_path)}", position=0
    ):
        if str(archive_path) in error_log["complete"]:
            continue

        start_time = time.time()

        try:
            if str(archive_path) in error_log["incomplete"]:
                print(f"Retrying {archive_path.name}...")

            merged_data = reader.process_file(
                archive_path,
                kmer_size=conf["db"]["kmer_size"],
                full=conf["db"]["full"],
                window_size=conf["db"]["window_size"],
                num_windows=conf["db"]["num_windows"],
            )
            writer.save_data(merged_data)

            error_log["complete"].append(str(archive_path))
            if str(archive_path) in error_log["incomplete"]:
                error_log["incomplete"].remove(str(archive_path))

            duration = time.time() - start_time
            error_log["duration"][archive_path.name] = format_duration(duration)

        except Exception as e:
            print(f"Error processing {archive_path.name}: {e}")
            traceback.print_exc()
            if str(archive_path) not in error_log["incomplete"]:
                error_log["incomplete"].append(str(archive_path))

        finally:
            with open(json_log_path, "w", encoding="utf-8") as error_log_file:
                json.dump(error_log, error_log_file, indent=4)


def load_config(json_file: Path | str):
    with open(json_file, "r") as file:
        conf = json.load(file)
    return conf


def get_current_datetime_string() -> str:
    return datetime.now().strftime("%Y_%m_%d_%H_%M_%S")


def train_model(conf: dict, rank: str | None, save_path: str | None):
    if rank not in RANKS:
        raise ValueError(f"Invalid rank {rank}")
    api = API(
        api_cache_dir=Path(conf["api"]["cache_dir"]),
        email=conf["api"]["email"],
        can_download=conf["api"]["can_download"],
    )
    model = XGBoostModel(api=api, use_gpu=conf["model"]["gpu"])
    db = Database(path=Path(conf["db"]["output_dir"]), api=api)
    mat = db.make_dmatrix(
        rank,
        sample_limit_by_tax_id=None,  # TODO
        normalize=conf["model"]["normalize"],
    )
    model.train(mat, tax_id_2_name=conf["model"]["tax_id_2_name"], kfold=None)
    if save_path is None:
        save_path = (
            Path(conf["model"]["default_models_dir"]) / get_current_datetime_string()
        ).with_suffix(".xgb")
    model.save(save_path)


def populate_api_cache(conf: dict):
    API(
        api_cache_dir=Path(conf["api"]["cache_dir"]),
        email=conf["api"]["email"],
        can_download=True,
    ).populate_api_cache(
        metadata_csv_path=Path(conf["allthebacteria"]["metadata_dir"])
        / METADATA_FILENAME
    )


def export_api_cache(conf: dict):
    API(
        api_cache_dir=Path(conf["api"]["cache_dir"]),
        email="",
        can_download=False,
    ).export_db()


def import_apt_cache(conf: dict):
    API(
        api_cache_dir=Path(conf["api"]["cache_dir"]),
        email="",
        can_download=False,
    ).import_db()


def debug():
    """debugging, ignore it"""

    mat = Database.deserialize_dmatrix("wisp_allthebacteria/out/mat2.pkl")
    api = API("/data/microtaxo/apicache", "cyrille.leroux@irisa.fr", True)
    model = XGBoostModel(api=api, use_gpu=False)
    res = model.train_with_kfold(mat, tax_id_2_name=True)
    print(res)

    # db = Database("/data/microtaxo/db_full_4", api)
    # mat = db.make_dmatrix("phylum", sample_limit_by_tax_id=None, normalize="min_max")
    # db.serialize_dmatrix(mat, "wisp_allthebacteria/out/mat2.pkl")
    # print((mat.num_row(), mat.num_col()))

    # api.export()
    # print(db.get_tax_ids_by_rank("phylum"))
    # data = db._get_tax_id_data(222)

    # writer = Writer(path="/data/microtaxo/db_full_4")
    # with open("/data/microtaxo/merged_data.pkl", "rb") as file:
    #     merged_data = pickle.load(file)
    # writer.save_data(merged_data)


if __name__ == "__main__":
    # debug()

    parser = argparse.ArgumentParser(description="AllTheBacteria Database Scripts")
    parser.add_argument(
        "--json",
        type=str,
        help="Path to the JSON file containing script configuration",
        default=Path(__file__).resolve().parent / "config/genouest.json",
    )
    parser.add_argument("--create-db", action="store_true", help="Create a database")
    parser.add_argument("--train", action="store_true", help="Train a model")

    parser.add_argument(
        "--populate-api-cache",
        action="store_true",
        help="get most of needed data from Entry (but you should use import/export api-cache instead)",
    )
    parser.add_argument(
        "--import-api-cache",
        action="store_true",
        help="pickle api cache to out/cache_dump.pkl",
    )
    parser.add_argument(
        "--export-api-cache",
        action="store_true",
        help="unpickle api cache from out/cache_dump.pkl",
    )

    parser.add_argument(
        "--train-model",
        action="store_true",
        help="Train a new model (may need --save-model)",
    )

    parser.add_argument(
        "--rank", type=str, help="Rank: phylum, kingdom, etc.", default=None
    )
    parser.add_argument(
        "--save-model",
        type=str,
        help="Path to save the trained model (or check config for default location)",
    )
    parser.add_argument(
        "--evaluate-model-kfolds",
        type=int,
        help="Number of folds for k-fold cross-validation (check config for report location)",
    )
    parser.add_argument(
        "--load-model",
        type=str,
        help="Path to load an existing model",
    )
    parser.add_argument(
        "--evaluate-fa",
        type=str,
        help="Path to a FASTA file for evaluation (check config for report location)",
    )

    args = parser.parse_args()

    conf = load_config(Path(args.json))

    if args.populate_api_cache:
        populate_api_cache(conf)

    if args.export_api_cache:
        export_api_cache(conf)

    if args.import_api_cache:
        import_apt_cache(conf)

    if args.create_db:
        create_db(conf)

    if args.train_model:
        train_model(conf=conf, rank=args.rank, save_path=args.save_model)

    # train-model, save-model, evaluate-model-kfolds, load-model, evaluate-fa
