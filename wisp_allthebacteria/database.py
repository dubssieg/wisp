import logging
from pathlib import Path
import sys
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
        self.clean()

    def clean(self):
        """undo unfinished transactions"""
        md_db = self._get_db(db_type="md")
        if CURRENT_TRANSACTION in md_db:
            LOG.warning(
                f"cleaning database {self.get_db_path()} (last transaction failed)"
            )
            data = md_db[CURRENT_TRANSACTION]
            merged_data = data["merged_data"]
            for tax_id in tqdm(merged_data.keys(), desc="cleaning cache"):
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
            LOG.warning("database cleaned")

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
        LOG.debug(f"file processed: {file_path}, start transaction")
        # from utils import deserialize
        # data = deserialize("wisp_allthebacteria/out/data.pkl")

        # mark transaction
        md_db[CURRENT_TRANSACTION] = data
        last_valid_ids = dict()

        merged_data = data["merged_data"]
        # warning, tax_id is a str
        LOG.debug(f"add data to: {self.get_db_path()}")
        for tax_id, tdata in tqdm(merged_data.items(), desc=f"Pushing {archive}"):
            tax_id = self._parse_tax_id(tax_id)
            last_valid_id = self._get_last_valid_id(tax_id)
            counters = tdata["counters"]
            sources = tdata["sources"]
            counter_db = self._get_db(db_type="counter", tax_id=tax_id)
            source_db = self._get_db(db_type="source", tax_id=tax_id)
            for i, (source, counter) in tqdm(
                enumerate(zip(sources, counters)),
                desc="Adding counters and sources",
                total=len(counters),
                leave=False,
            ):
                current_id = last_valid_id + i + 1
                counter_db[current_id] = counter
                source_db[current_id] = source
                last_valid_ids[tax_id] = current_id

        LOG.debug(f"data added: {self.get_db_path()}, ending transaction")
        # end transaction
        del md_db[CURRENT_TRANSACTION]
        for tax_id, last_valid_id in last_valid_ids.items():
            self._set_last_valid_id(tax_id=tax_id, last_valid_id=last_valid_id)
        archives.append(archive)
        md_db[ARCHIVES] = archives
        LOG.debug("transaction ended successfully")

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
