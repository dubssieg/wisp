import logging
from pathlib import Path
import sys
import threading
from typing import Literal

from diskcache import Cache
from functools import lru_cache
from tqdm.auto import tqdm

from reader import Reader

LOG = logging.getLogger(__name__)

DB_TYPE = Literal["counter", "source", "index", "md"]
LAST_VALID_ID = "_last_valid_id_"
CURRENT_TRANSACTION = "_current_transaction_"
ARCHIVES = "_archives_"


class Database:
    def __init__(
        self,
        kmer_size: int,
        window_size: int,
        step: int,
        full: bool,
        dbs_path: str | Path,
        reader: Reader,
    ):
        self._dbs_path = Path(dbs_path).resolve()
        self._kmer_size = kmer_size
        self._reader = reader
        self._window_size = window_size
        self._step = step
        self._full = full
        self._db_lock = threading.Lock()
        self.clean()

    def get_info(self) -> dict:
        """Get DB infos: tax_ids, number of counters, etc."""
        # tax_ids
        counters_dir = self.get_db_path() / "counter"
        tax_ids = [
            int(dir.name)
            for dir in counters_dir.iterdir()
            if dir.is_dir() and dir.name.isdigit()
        ]

        # archives pushed
        md_db = self._get_db(db_type="md")
        archives = md_db.get(ARCHIVES, [])

        # counters for each tax_id
        counters = {}
        for tax_id in tax_ids:
            counter_db = self._get_db(db_type="counter", tax_id=tax_id)
            counters[tax_id] = len(counter_db)

        return {"tax_ids": tax_ids, "archives": archives, "counters": counters}

    def clean(self):
        """Undo unfinished transactions"""
        md_db = self._get_db(db_type="md")
        if CURRENT_TRANSACTION in md_db:
            LOG.warning(
                f"Cleaning database {self.get_db_path()} (last transaction failed)"
            )
            data = md_db[CURRENT_TRANSACTION]
            merged_data = data["merged_data"]
            for tax_id in tqdm(merged_data.keys(), desc="cleaning DB"):
                counter_db = self._get_db(db_type="counter", tax_id=tax_id)
                source_db = self._get_db(db_type="source", tax_id=tax_id)
                current_id = self._get_last_valid_id(tax_id)
                deleting = True
                while deleting:
                    current_id += 1
                    deleting = False
                    if current_id in counter_db:
                        del counter_db[current_id]
                        deleting = True
                    if current_id in source_db:
                        del source_db[current_id]
                        deleting = True
            del md_db[CURRENT_TRANSACTION]
            LOG.warning("Database cleaned")

    def push_file(self, file_path: str | Path):
        """Add content."""
        file_path = Path(file_path).resolve()
        LOG.debug(f"push file: {file_path}")
        md_db = self._get_db(db_type="md")
        archive = file_path.stem.split(".")[0]
        archives = md_db.get(ARCHIVES, [])

        if archive in archives:
            LOG.info(f"{archive} already in DB: skip")
            return

        data = self._reader.process_file(
            file_path=file_path,
            kmer_size=self._kmer_size,
            window_size=self._window_size,
            step=self._step,
            full=self._full,
        )
        LOG.debug(f"File processed: {file_path}, start transaction")

        db_thread = threading.Thread(target=self._add_data_to_db, kwargs={"data": data})
        db_thread.start()

    def _add_data_to_db(self, data):
        with self._db_lock:
            md_db = self._get_db(db_type="md")
            md_db[CURRENT_TRANSACTION] = data
            archives = md_db.get(ARCHIVES, [])
            archive = data["archive"].stem.split(".")[0]

            last_valid_ids = {}

            merged_data = data["merged_data"]
            # warning, tax_id is a str
            LOG.debug(f"Add data to: {self.get_db_path()}")
            for tax_id, tdata in tqdm(
                merged_data.items(), desc=f"Pushing {archive}", position=1, leave=False
            ):
                tax_id = self._parse_tax_id(tax_id)
                last_valid_id = self._get_last_valid_id(tax_id)
                counters = tdata["counters"]
                sources = tdata["sources"]
                counter_db = self._get_db(db_type="counter", tax_id=tax_id)
                source_db = self._get_db(db_type="source", tax_id=tax_id)

                # populate and add batch of counters/sources
                batch_counters = {}
                batch_sources = {}

                for i, (source, counter) in tqdm(
                    enumerate(zip(sources, counters)),
                    desc="Adding counters and sources",
                    total=len(counters),
                    position=2,
                    leave=False,
                ):
                    current_id = last_valid_id + i + 1
                    batch_counters[current_id] = counter
                    batch_sources[current_id] = source

                LOG.debug(f"Pushing {len(batch_counters)} counters to DB")
                with counter_db.transact():
                    for current_id, counter in tqdm(
                        batch_counters.items(), desc="Counters", position=3, leave=False
                    ):
                        counter_db[current_id] = counter
                LOG.debug(f"Pushing {len(batch_sources)} sources to DB")
                with source_db.transact():
                    for current_id, source in tqdm(
                        batch_sources.items(), desc="Sources", position=3, leave=False
                    ):
                        source_db[current_id] = source
                LOG.debug("Counters and sources successfully added")

            LOG.debug(f"Data added: {self.get_db_path()}, ending transaction")
            # end transaction
            del md_db[CURRENT_TRANSACTION]
            for tax_id, last_valid_id in last_valid_ids.items():
                self._set_last_valid_id(tax_id=tax_id, last_valid_id=last_valid_id)
            archives.append(archive)
            md_db[ARCHIVES] = archives
            LOG.debug("Transaction ended successfully")

    def _get_last_valid_id(self, tax_id: int) -> int:
        """Get last inserted id"""
        md_db = self._get_db(db_type="md", tax_id=tax_id)
        return md_db.get(LAST_VALID_ID, -1)

    def _set_last_valid_id(self, tax_id: int, last_valid_id: int) -> int:
        """Get last inserted id"""
        md_db = self._get_db(db_type="md", tax_id=tax_id)
        last_valid_id = md_db.get(LAST_VALID_ID)
        if last_valid_id is None:
            return -1
        return last_valid_id

    def _parse_tax_id(self, tax_id: str | int) -> str | int:
        try:
            return int(tax_id)
        except (ValueError, TypeError):
            return tax_id

    @lru_cache(maxsize=100)
    def _get_db(self, db_type: DB_TYPE, tax_id: int | None = None) -> Cache:
        """Get sub DB"""
        path = self.get_db_path()
        path /= db_type
        if tax_id:
            path /= str(tax_id)

        path.mkdir(parents=True, exist_ok=True)

        return Cache(path, size_limit=sys.maxsize)

    def get_db_path(self):
        path = self._dbs_path / str(self._kmer_size)
        if self._full:
            path /= "full"
        else:
            path /= f"{self._window_size}_{self._step}"
        return path
