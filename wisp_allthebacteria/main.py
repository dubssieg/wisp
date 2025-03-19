import gc
import json
import logging
import os
import sys
import time
from pathlib import Path
import argparse
from tqdm.auto import tqdm
from metadata import Metadata
from reader import Reader
from api import API
from model import XGBoostModel
from database import Database, DatabaseBuilder
from utils import (
    SystemStatsLogger,
    format_duration,
    get_current_datetime_string,
    cpu_count,
    config_logger,
)

LOG = logging.getLogger(__name__)


METADATA_FILENAME = "ena_metadata.tsv"

RANKS = [
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "subclass",
    "order",
    "suborder",
    "family",
    "subfamily",
    "tribe",
    "genus",
    "subgenus",
    "species",
    "subspecies",
    "strain",
    # ?
    "biotype",
    "clade",
    "no rank",
    "pathogroup",
    "serogroup",
    "serotype",
    "species group",
    "species subgroup",
]


def create_db(conf: dict):
    LOG.info("create_db")
    with SystemStatsLogger(interval=30, pid=os.getpid()):  # TODO: conf + level
        input_path = Path(conf["allthebacteria"]["assembly_dir"])
        output_path = Path(conf["db"]["path"])
        metadata_path = Path(conf["allthebacteria"]["metadata_dir"]) / METADATA_FILENAME
        api_cache_path = Path(conf["api"]["cache_dir"])
        output_path.mkdir(parents=True, exist_ok=True)

        archives = list(input_path.glob("*.xz"))

        api = API(
            api_cache_dir=api_cache_path,
            email=conf["api"]["email"],
            can_download=conf["api"]["can_download"],
            preload=False,
        )
        md = Metadata(csv_path=metadata_path, api=api, start_loaded=True)
        num_workers = cpu_count(conf["db"]["create_db_workers"])
        num_threads = cpu_count(conf["db"]["db_insert_threads"])
        LOG.info(
            f"Max CPUs: {cpu_count('max')}, using {num_workers} fasta workers and {num_threads} DB insertion threads"
        )
        reader = Reader(md, num_workers=num_workers)

        db = DatabaseBuilder(
            kmer_size=conf["db"]["kmer_size"],
            window_size=conf["db"]["window_size"],
            step=conf["db"]["step"],
            full=conf["db"]["full"],
            dbs_path=output_path,
            fanout_shards=conf["db"]["fanout_shards"],
            reader=reader,
            fasta_batch_size=conf["db"]["fasta_batch_size"],
            insert_threads=num_threads,
            compressed=conf["db"]["compressed"],
            merged_data_as_db=conf["db"]["merged_data_as_db"],
        )

        json_create_db_path = db.get_db_path() / "create_db.json"
        json_create_db_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            json_create_db = json.loads(json_create_db_path.read_text(encoding="utf-8"))
        except (FileNotFoundError, json.JSONDecodeError):
            json_create_db = {"complete": [], "incomplete": [], "duration": {}}

        complete = set(json_create_db["complete"])
        incomplete = set(json_create_db["incomplete"])
        duration = json_create_db["duration"]

        try:
            for i, archive_path in tqdm(
                enumerate(archives),
                desc=f"Processing {input_path} -> {output_path}",
                total=len(archives),
            ):
                LOG.debug(f"=== [{i+1} / {len(archives)}] {archive_path.name} ===")
                archive_str = str(archive_path)

                start_time = time.time()

                try:
                    if archive_str in incomplete:
                        LOG.info(f"Retrying {archive_path.name}...")

                    if db.has_archive(archive_path):
                        LOG.warning(f"{archive_path.name} already in DB, skipping.")
                    else:

                        db.push_file(archive_path)
                        gc.collect()

                    complete.add(archive_str)
                    incomplete.discard(archive_str)

                    duration[archive_path.name] = format_duration(
                        time.time() - start_time
                    )

                except Exception:
                    LOG.exception(f"Error processing {archive_path.name}")
                    incomplete.add(archive_str)
                    raise

                finally:
                    with json_create_db_path.open("w", encoding="utf-8") as json_file:
                        json.dump(
                            {
                                "complete": list(complete),
                                "incomplete": list(incomplete),
                                "duration": duration,
                            },
                            json_file,
                            indent=4,
                        )
        except Exception:
            LOG.exception(f"Error processing {archive_path.name}")
            raise

    LOG.info("create_db -> DONE")


def db_info(conf: dict):
    LOG.info("DB info")

    db = Database(
        kmer_size=conf["db"]["kmer_size"],
        window_size=conf["db"]["window_size"],
        step=conf["db"]["step"],
        full=conf["db"]["full"],
        dbs_path=conf["db"]["path"],
        compressed=conf["db"]["compressed"],
        fanout_shards=conf["db"]["fanout_shards"],
    )

    print(db.get_info(as_str=True))


def load_config(json_file: Path | str):
    with open(json_file, "r") as file:
        conf = json.load(file)
    return conf


def train_model(conf: dict, rank: str | None, save_path: str | None, kfold: int | None):
    LOG.info("train_model")
    if rank not in RANKS:
        raise ValueError(f"Invalid rank {rank}")

    if kfold is not None and kfold == -1:
        kfold = conf["model"]["kfold"]

    api = API(
        api_cache_dir=Path(conf["api"]["cache_dir"]),
        email=conf["api"]["email"],
        can_download=conf["api"]["can_download"],
    )
    model = XGBoostModel(
        api=api,
        use_gpu=conf["model"]["gpu"],
        scientific_name=conf["model"]["scientific_name"],
    )
    db = Database(path=Path(conf["db"]["output_dir"]), api=api)
    batch_size = conf["model"]["batch_size"]
    # dmat_generator = db.make_dmatrix(
    #     rank=rank,
    #     sample_limit_by_tax_id=None,  # TODO
    #     normalize=conf["model"]["normalize"],
    #     batch_size=batch_size,
    # )
    # dgf = DMatrixGeneratorFactory(
    #     database=db,
    #     rank=rank,
    #     normalize=conf["model"]["normalize"],
    #     sample_limit_by_tax_id=None,
    #     batch_size=batch_size,
    #     max_samples=conf["model"]["max_samples"],
    # )
    dmat = db.make_dmatrix(
        rank=rank,
        normalize=conf["model"]["normalize"],
    )
    tax_id_classes = db.get_tax_id_classes(rank)

    report = model.train(
        # dtrain=dmat_generator,
        # dtrain_factory=dgf,
        dtrain=dmat,
        tax_id_classes=tax_id_classes,
        kfold=kfold,
        num_boost_round=conf["model"]["num_boost_round"],
    )
    # train mode
    if kfold is None:
        if save_path is None:
            save_path = (
                Path(conf["model"]["default_models_dir"])
                / get_current_datetime_string()
            )
        model.save(save_path)
    # evaluate mode with kfold
    else:
        if save_path is None:
            save_path = (
                Path(conf["report"]["default_reports_dir"])
                / get_current_datetime_string()
            )
    report_header = {"header": {"Rang": rank, "batch_size": batch_size}}

    model.save_report(dir_path=save_path, additional_data=report_header)
    print(report)


# def evaluate_model(conf: dict, nfolds: int, model_path: str|Path|None):
#     if model_path is None:
#             raise ValueError("Need --load-model argument")

#     api = API(
#         api_cache_dir=Path(conf["api"]["cache_dir"]),
#         email=conf["api"]["email"],
#         can_download=conf["api"]["can_download"],
#     )
#     model = XGBoostModel(api=api, use_gpu=conf["model"]["gpu"])
#     model.load(model_path)


def populate_api_cache(conf: dict):
    LOG.info("populate_api_cache")
    API(
        api_cache_dir=Path(conf["api"]["cache_dir"]),
        email=conf["api"]["email"],
        can_download=True,
    ).populate_api_cache(
        metadata_csv_path=Path(conf["allthebacteria"]["metadata_dir"])
        / METADATA_FILENAME
    )


def export_api_cache(conf: dict):
    LOG.info("export_apt_cache")
    API(
        api_cache_dir=Path(conf["api"]["cache_dir"]),
        email="",
        can_download=False,
    ).export_db()


def import_apt_cache(conf: dict):
    LOG.info("import_apt_cache")
    API(
        api_cache_dir=Path(conf["api"]["cache_dir"]),
        email="",
        can_download=False,
    ).import_db()


def debug(conf):
    """debugging, ignore it"""
    LOG.info("debug")

    api = API(
        api_cache_dir=conf["api"]["cache_dir"],
        email=conf["api"]["email"],
        can_download=conf["api"]["can_download"],
        preload=False,
    )

    db = Database(
        kmer_size=conf["db"]["kmer_size"],
        window_size=conf["db"]["window_size"],
        step=conf["db"]["step"],
        full=conf["db"]["full"],
        dbs_path=conf["db"]["path"],
        fanout_shards=conf["db"]["fanout_shards"],
        compressed=conf["db"]["compressed"],
        api=api,
    )

    # print(db.get_info(as_str=True))
    tids = db.get_tax_ids_by_rank("phylum")
    res = dict()
    for rank_tax_id, tax_ids in tids.items():
        sids = db.tax_ids_to_sample_ids(tax_ids)
        res[rank_tax_id] = sids

    print({k: len(v) for k, v in res.items()})
    pass

    # mat = Database.deserialize_dmatrix("wisp_allthebacteria/out/mat2.pkl")

    # api = API("/data/microtaxo/apicache", "cyrille.leroux@irisa.fr", True)
    # db = Database("/data/microtaxo/db_full_4", api)
    # db.index_by_rank("phylum")
    # matgen = db.make_dmatrix(rank="phylum", normalize="min_max", batch_size=100)
    # tax_id_classes = db.get_tax_id_classes("phylum")
    # model = XGBoostModel(api=api, use_gpu=False)
    # res = model.train(matgen, kfold=None, tax_id_classes=tax_id_classes)

    # report_header = {
    #     "header": {
    #         "Rang": "phylum",
    #     }
    # }
    # model.save_report(
    #     dir_path="wisp_allthebacteria/out/report1", additional_data=report_header
    # )
    # print(res)

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

    parser = argparse.ArgumentParser(
        description="AllTheBacteria Database Scripts",
        epilog="""
        Examples:
        Create a database:
            python main.py --create-db
        
        Create a database (local machine/debug)
            python main.py -- create-db --json="wisp_allthebacteria/config/clx_debug.json"

        Train a model for phylum classification
            python main.py --train-model --rank="phylum"

        Evaluate a model for phylum classification 1
            python main.py --evaluate-model-kfold --rank="phylum"
        
        Evaluate a model for phylum classification 2
            python main.py --evaluate-model-kfold=5 --rank="phylum" --save-path="report_1"

        """,
    )
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
        "--save-path",
        type=str,
        help="Path to save results (or check config for default location)",
    )

    parser.add_argument(
        "--evaluate-model-kfolds",
        type=int,
        nargs="?",
        const=-1,
        help="Number of folds for k-fold cross-validation (check config for report location)",
    )

    parser.add_argument(
        "--load-model", type=str, help="Path to load an existing model", default=None
    )
    parser.add_argument(
        "--evaluate-fa",
        type=str,
        help="Path to a FASTA file for evaluation (check config for report location)",
    )

    parser.add_argument(
        "--db-info",
        action="store_true",
        help="Show database informations",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Dev only, do not use",
    )

    args = parser.parse_args()

    conf = load_config(Path(args.json))

    config_logger(**conf["log"])

    if args.debug:
        debug(conf)
        sys.exit("debug")

    if args.populate_api_cache:
        populate_api_cache(conf)

    if args.export_api_cache:
        export_api_cache(conf)

    if args.import_api_cache:
        import_apt_cache(conf)

    if args.create_db:
        create_db(conf)

    if args.db_info:
        db_info(conf)

    if args.train_model:
        train_model(conf=conf, rank=args.rank, save_path=args.save_path, kfold=None)

    if args.evaluate_model_kfolds:
        train_model(
            conf=conf,
            rank=args.rank,
            save_path=args.save_path,
            kfold=args.evaluate_model_kfolds,
        )

    # train-model, save-model, evaluate-model-kfolds, load-model, evaluate-fa
