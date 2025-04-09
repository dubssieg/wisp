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
from fastakmer import FastaKmer
from taxdb import TaxDB, RANKS
from model import XGBoostModel
from supermodel import SuperModel
from database import Database, DatabaseBuilder
from dataset import Dataset
from utils import (
    SystemStatsLogger,
    format_duration,
    get_current_datetime_string,
    cpu_count,
    config_logger,
)

LOG = logging.getLogger(__name__)


METADATA_FILENAME = "ena_metadata.tsv"


def create_db(conf: dict):
    LOG.info("Create DB")
    with SystemStatsLogger(interval=30, pid=os.getpid()):  # TODO: conf + level
        input_path = Path(conf["allthebacteria"]["assembly_dir"])
        output_path = Path(conf["db"]["path"])
        metadata_path = Path(conf["allthebacteria"]["metadata_dir"]) / METADATA_FILENAME
        output_path.mkdir(parents=True, exist_ok=True)

        taxdb = TaxDB(
            cache_dir=conf["taxdb"]["cache_dir"],
            email=conf["taxdb"]["email"],
            can_download=conf["taxdb"]["can_download"],
            preload=False,
        )
        md = Metadata(csv_path=metadata_path, taxdb=taxdb, start_loaded=True)
        num_workers = cpu_count(conf["db"]["create_db_workers"])
        num_threads = cpu_count(conf["db"]["db_insert_threads"])
        LOG.info(
            f"Max CPUs: {cpu_count('max')}, using {num_workers} fasta workers and {num_threads} DB insertion threads"
        )

        fasta_kmer = FastaKmer(
            md,
            num_workers=num_workers,
            taxdb=taxdb,
            tmp_path=conf["db"]["tmp_path"],
        )
        db = DatabaseBuilder(
            kmer_sizes=conf["db"]["kmer_sizes"],
            window_size=conf["db"]["window_size"],
            step=conf["db"]["step"],
            full=conf["db"]["full"],
            dbs_path=output_path,
            fanout_shards=conf["db"]["fanout_shards"],
            fasta_kmer=fasta_kmer,
            fasta_batch_size=conf["db"]["fasta_batch_size"],
            insert_threads=num_threads,
            compression=conf["db"]["compression"],
            merged_data_as_db=conf["db"]["merged_data_as_db"],
            species_count_limit=conf["db"]["species_count_limit"],
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

        archives = [
            archive
            for archive in input_path.glob("*.xz")
            if not db.has_archive(archive)
        ]

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

                    db.push_archive(archive_path)
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
        kmer_sizes=conf["db"]["kmer_sizes"],
        window_size=conf["db"]["window_size"],
        step=conf["db"]["step"],
        full=conf["db"]["full"],
        dbs_path=conf["db"]["path"],
        compression=conf["db"]["compression"],
        fanout_shards=conf["db"]["fanout_shards"],
    )

    print(db.get_info(as_str=True))


def batch_count(conf: dict, rank: str):
    LOG.info("Batch counts")
    database = Database(
        kmer_sizes=conf["db"]["kmer_sizes"],
        window_size=conf["db"]["window_size"],
        step=conf["db"]["step"],
        full=conf["db"]["full"],
        dbs_path=conf["db"]["path"],
        compression=conf["db"]["compression"],
        fanout_shards=conf["db"]["fanout_shards"],
    )
    taxdb = TaxDB(
        cache_dir=conf["taxdb"]["cache_dir"],
        email=conf["taxdb"]["email"],
        can_download=conf["taxdb"]["can_download"],
        preload=False,
    )
    ds = Dataset(database=database, taxdb=taxdb)
    batch_size = conf["model"]["batch_size"]
    print(f"Batch size= {batch_size}\n")
    if min_samples := conf["model"]["min_samples_per_class"]:

        print(f"WITH Min samples by class = {min_samples}\n")
        print(
            ds.get_batch_counts_for_rank(
                rank=rank,
                batch_size=batch_size,
                min_samples_per_class=conf["model"]["min_samples_per_class"],
                as_str=True,
            )
        )

    print("\n\nWITHOUT Min samples by class\n")
    print(
        ds.get_batch_counts_for_rank(
            rank=rank,
            batch_size=batch_size,
            min_samples_per_class=None,
            as_str=True,
        )
    )


def get_counters(conf: dict, tax_id: int, limit: int | None):

    taxdb = TaxDB(
        cache_dir=conf["taxdb"]["cache_dir"],
        email=conf["taxdb"]["email"],
        can_download=conf["taxdb"]["can_download"],
        preload=False,
    )
    name = taxdb[tax_id]["ScientificName"]
    name = name.replace(" ", "_").lower()
    path = Path(name).with_suffix(".csv")

    database = Database(
        kmer_sizes=conf["db"]["kmer_sizes"],
        window_size=conf["db"]["window_size"],
        step=conf["db"]["step"],
        full=conf["db"]["full"],
        dbs_path=conf["db"]["path"],
        compression=conf["db"]["compression"],
        fanout_shards=conf["db"]["fanout_shards"],
    )
    df = database.get_tax_id_samples(tax_id, as_dataframe=True, limit=limit)

    LOG.info(f"Writing file {path}...")
    df.to_csv(path, index=False)


def load_config(json_file: Path | str):
    with open(json_file, "r") as file:
        conf = json.load(file)
    return conf


def train_model(
    conf: dict,
    rank: str | None,
    save_path: str | None,
    kfold: bool,
):
    LOG.info("Train")
    with SystemStatsLogger(interval=30):
        if rank not in RANKS:
            raise ValueError(f"Invalid rank {rank}")

        taxdb = TaxDB(
            cache_dir=conf["taxdb"]["cache_dir"],
            email=conf["taxdb"]["email"],
            can_download=conf["taxdb"]["can_download"],
        )
        database = Database(
            kmer_sizes=conf["db"]["kmer_sizes"],
            window_size=conf["db"]["window_size"],
            step=conf["db"]["step"],
            full=conf["db"]["full"],
            dbs_path=conf["db"]["path"],
            fanout_shards=conf["db"]["fanout_shards"],
            compression=conf["db"]["compression"],
        )

        normalize = conf["model"]["normalize"]
        if not normalize:
            normalize = None
        dt = get_current_datetime_string()
        if conf["model"]["force_workspaces_dir"]:
            workspace_path = Path(conf["model"]["force_workspaces_dir"]).resolve()
        else:
            workspace_path = (
                Path(conf["model"]["default_workspaces_dir"]).resolve() / dt
            )

        xgb_model = XGBoostModel(
            rank=rank,
            database=database,
            normalize=normalize,
            batch_size=conf["model"]["batch_size"],
            taxdb=taxdb,
            generator_threads=cpu_count(conf["model"]["generator_threads"]),
            sample_balance_factor=conf["model"]["sample_balance_factor"],
            batch_balance_factor=conf["model"]["batch_balance_factor"],
            min_samples_per_class=conf["model"]["min_samples_per_class"],
            max_buffer_total_size=conf["model"]["max_buffer_total_size"],
            workspace_path=workspace_path,
            start_generator=True,
        )

        train_batch_count = conf["model"]["train_batch_count"]
        if train_batch_count == "max":
            train_batch_count = None
        if save_path is None:
            save_path = Path(conf["model"]["default_models_dir"]) / dt

        if kfold:
            xgb_model.kfold(
                k=conf["model"]["kfold"],
                save_path=save_path,
                train_batch_count=train_batch_count,
                eval_patience=conf["model"]["eval_patience"],
                eval_batch_count=conf["model"]["eval_batch_count"],
                num_boost_round=conf["model"]["num_boost_round"],
            )
        else:
            xgb_model.train(
                save_path=save_path,
                train_batch_count=train_batch_count,
                eval_patience=conf["model"]["eval_patience"],
                eval_batch_count=conf["model"]["eval_batch_count"],
                test_batch_count=conf["model"]["test_batch_count"],
                num_boost_round=conf["model"]["num_boost_round"],
            )

        xgb_model.stop()


def train_supermodel(conf, name: str):
    LOG.info("Train supermodel")

    taxdb = TaxDB(
        cache_dir=conf["taxdb"]["cache_dir"],
        email=conf["taxdb"]["email"],
        can_download=conf["taxdb"]["can_download"],
    )

    database = Database(
        kmer_sizes=conf["db"]["kmer_sizes"],
        window_size=conf["db"]["window_size"],
        step=conf["db"]["step"],
        full=conf["db"]["full"],
        dbs_path=conf["db"]["path"],
        fanout_shards=conf["db"]["fanout_shards"],
        compression=conf["db"]["compression"],
    )

    sm = SuperModel(
        model_conf=conf["model"],
        supermodel_conf=conf["supermodels"][name],
        database=database,
        taxdb=taxdb,
    )
    sm.train()


def populate_taxdb(conf: dict):
    LOG.info("Populate_taxdb")
    TaxDB(
        cache_dir=Path(conf["taxdb"]["cache_dir"]),
        email=conf["taxdb"]["email"],
        can_download=True,
    ).populate_api_cache(
        metadata_csv_path=Path(conf["allthebacteria"]["metadata_dir"])
        / METADATA_FILENAME
    )


def export_taxdb(conf: dict):
    LOG.info("Export_taxdb")
    TaxDB(
        cache_dir=Path(conf["taxdb"]["cache_dir"]),
        email="",
        can_download=False,
    ).export_db()


def import_taxdb(conf: dict):
    LOG.info("Import_taxdb")
    TaxDB(
        cache_dir=Path(conf["taxdb"]["cache_dir"]),
        email="",
        can_download=False,
    ).import_db()


def debug(conf):
    """debugging, ignore it"""
    LOG.info("Debug")

    taxdb = TaxDB(
        cache_dir=conf["taxdb"]["cache_dir"],
        email=conf["taxdb"]["email"],
        can_download=conf["taxdb"]["can_download"],
    )

    database = Database(
        kmer_sizes=conf["db"]["kmer_sizes"],
        window_size=conf["db"]["window_size"],
        step=conf["db"]["step"],
        full=conf["db"]["full"],
        dbs_path=conf["db"]["path"],
        fanout_shards=conf["db"]["fanout_shards"],
        compression=conf["db"]["compression"],
    )

    # dataset = Dataset(database=database, taxdb=taxdb)

    sm = SuperModel(
        model_conf=conf["model"],
        supermodel_conf=conf["supermodels"]["model1"],
        database=database,
        taxdb=taxdb,
    )
    sm.train()

    # rank_an = ds.get_ranks_analysis(scientific_names=False)
    # res = {}
    # for rank_tid in rank_an["kingdom"]["tax_ids"].keys():
    #     res[rank_tid] = ds._get_tax_ids_by_rank(
    #         "phylum", parent_filter={"kingdom": rank_tid}
    #     )

    # print(res)
    # print(ds._get_tax_ids_by_rank("phylum"))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="AllTheBacteria Database Scripts",
        epilog="""
        Examples:
        Create a database:
            python main.py --create-db

        Create a database (local machine/debug)
            python main.py -- create-db --conf="config/clx_debug.json"

        Evaluate a model for phylum classification
            python main.py --train --kfold --rank="phylum"

        Train a model for phylum classification
            python main.py --train --rank="phylum"
        """,
    )
    parser.add_argument(
        "--conf",
        type=str,
        help="Path to the JSON file containing script configuration",
        default=Path(__file__).resolve().parent / "config/genouest.json",
    )
    parser.add_argument("--create-db", action="store_true", help="Create a database")
    parser.add_argument("--train-model", action="store_true", help="Train a model")
    parser.add_argument("--train-supermodel", type=str, help="Train a supermodel")

    parser.add_argument(
        "--populate-taxdb",
        action="store_true",
        help="get most of needed data from Entrez",
    )
    parser.add_argument(
        "--import-taxdb",
        action="store_true",
        help="pickle taxdb to out/cache_dump.pkl",
    )
    parser.add_argument(
        "--export-taxdb",
        action="store_true",
        help="unpickle taxdb from out/cache_dump.pkl",
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
        "--kfold",
        action="store_true",
        help="Number of folds for k-fold cross-validation",
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
        "--batch-count",
        type=str,
        help="How many batches for this rank (use batch_size and min_samples_per_class). Ex: --batch_counts=phylum",
    )

    parser.add_argument(
        "--get-counters",
        type=int,
        help="Get counters for this tax-id. Ex: --get-counters=1464 --limit=30000",
        default=None,
    )

    parser.add_argument(
        "--limit",
        type=int,
        help="Limit (--get-counters only for now)",
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="Dev only, do not use",
    )

    args = parser.parse_args()

    conf = load_config(Path(args.conf))

    config_logger(**conf["log"])

    if args.debug:
        debug(conf)
        sys.exit("debug")

    if args.populate_taxdb:
        populate_taxdb(conf)

    if args.export_taxdb:
        export_taxdb(conf)

    if args.import_taxdb:
        import_taxdb(conf)

    if args.create_db:
        create_db(conf)

    if args.db_info:
        db_info(conf)

    if args.batch_count:
        batch_count(conf, rank=args.batch_count)

    if args.get_counters:
        get_counters(conf=conf, tax_id=args.get_counters, limit=args.limit)

    if args.train_model:
        train_model(
            conf=conf,
            rank=args.rank,
            save_path=args.save_path,
            kfold=args.kfold,
        )
    if args.train_supermodel:
        train_supermodel(conf=conf, name=args.train_supermodel)

    # train-model, save-model, evaluate-model-kfolds, load-model, evaluate-fa
