import json
import time
import traceback
from pathlib import Path

from tqdm.auto import tqdm
from metadata import Metadata
from reader import Reader
from writer import Writer
from api import API
import argparse


METADATA_FILENAME = "ena_metadata.tsv"


def format_duration(seconds):
    """Convert seconds to a formatted duration string (hh:mm:ss)."""
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{int(hours):02}:{int(minutes):02}:{int(seconds):02}"


def create_db(
    input_path: Path | str,
    output_path: Path | str,
    metadata_path: Path | str,
    api_cache_path: Path | str,
    api_can_download: bool,
    email: str,
    kmer_size: int,
    full: int,
    window_size: int,
    num_windows: int,
):
    input_path = Path(input_path)
    output_path = Path(output_path)
    metadata_path = Path(metadata_path)
    api_cache_path = Path(api_cache_path)
    output_path.mkdir(parents=True, exist_ok=True)

    json_log_path = output_path / "error_log.json"

    if json_log_path.exists():
        with open(json_log_path, "r", encoding="utf-8") as error_log_file:
            error_log = json.load(error_log_file)
    else:
        error_log = {"complete": [], "incomplete": [], "duration": {}}

    archives = list(input_path.glob("*.xz"))

    api = API(api_cache_dir=api_cache_path, email=email, can_download=api_can_download)
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
                kmer_size=kmer_size,
                full=full,
                window_size=window_size,
                num_windows=num_windows,
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


def debug():
    """debugging, ignore it"""
    from database import Database

    mat1 = Database.deserialize_dmatrix("wisp_allthebacteria/out/mat.pkl")

    api = API("/data/microtaxo/apicache", "cyrille.leroux@irisa.fr", True)
    db = Database("/data/microtaxo/db_light_4", api)
    mat = db.make_dmatrix("phylum", sample_limit_by_tax_id=None, normalize="min_max")
    db.serialize_dmatrix(mat, "wisp_allthebacteria/out/mat.pkl")
    print(mat.shape)

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

    args = parser.parse_args()

    conf = load_config(Path(args.json))

    if args.export_api_cache:
        API(
            api_cache_dir=Path(conf["api"]["cache_dir"]),
            email=conf["api"]["email"],
            can_download=True,
        ).populate_api_cache(
            metadata_csv_path=Path(conf["allthebacteria"]["metadata_dir"])
            / METADATA_FILENAME
        )

    if args.export_api_cache:
        API(
            api_cache_dir=Path(conf["api"]["cache_dir"]),
            email="",
            can_download=False,
        ).export_db()

    if args.import_api_cache:
        API(
            api_cache_dir=Path(conf["api"]["cache_dir"]),
            email="",
            can_download=False,
        ).import_db()

    if args.create_db:
        create_db(
            input_path=Path(conf["allthebacteria"]["assembly_dir"]),
            output_path=Path(conf["db"]["output_path"]),
            metadata_path=Path(conf["allthebacteria"]["metadata_dir"])
            / METADATA_FILENAME,
            api_cache_path=Path(conf["api"]["cache_dir"]),
            api_can_download=conf["api"]["can_download"],
            email=conf["api"]["email"],
            kmer_size=conf["db"]["kmer_size"],
            full=conf["db"]["full"],
            window_size=conf["db"]["window_size"],
            num_windows=conf["db"]["num_windows"],
        )
