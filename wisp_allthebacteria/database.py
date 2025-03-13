import logging
from pathlib import Path
import sys
import threading
from typing import Literal

from diskcache import Cache
from functools import lru_cache

from reader import Reader
from utils import slurm_tqdm, space_format

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

    def get_info(self, as_str=False) -> dict:
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

        # path
        info = {"tax_ids": tax_ids, "archives": archives, "counters": counters}

        if as_str:
            info_lines = []
            info_lines.append("\n=== Paths/Config ===")
            info_lines.append(f"Current base path: {self.get_db_path()}")
            info_lines.append("Available sub DBs:")
            all_dbs = self.list_dbs(self._dbs_path)
            for kmer, win_step in all_dbs:
                win, step = win_step.split("_")
                info_lines.append(
                    f"  - k-mer size: {kmer}, window size: {win}, step: {step}"
                )

            info_lines.append("\n=== Sequences ===")
            info_lines.append(f"{len(tax_ids)} counters -> [tax_id] sample_count:")
            items_per_line = 5
            counters_str = {
                tax_id: space_format(num) for tax_id, num in sorted(counters.items())
            }
            sorted_items = sorted(counters_str.items())
            max_tax_id_len = max(len(str(tax_id)) for tax_id, _ in sorted_items)
            max_counter_len = max(len(str(counter)) for _, counter in sorted_items)
            for i in range(0, len(sorted_items), items_per_line):
                line_items = sorted_items[i : i + items_per_line]
                line = "  - " + " | ".join(
                    [
                        f"[{tax_id:<{max_tax_id_len}}] {counter:<{max_counter_len}}"
                        for tax_id, counter in line_items
                    ]
                )
                info_lines.append(line)

            info_lines.append("")
            info_lines.append(f"Total: {space_format(sum(counters.values()))} samples")

            info_lines.append("\n==== Archives pushed ====")
            info_lines.append(", ".join(sorted(archives)))

            return "\n".join(info_lines)

        return info

    def has_archive(self, path: str | Path) -> bool:
        """Check if archive already was pushed in base."""
        md_db = self._get_db(db_type="md")
        return self._archive_stem(path) in md_db.get(ARCHIVES, [])

    def clean(self):
        """Undo unfinished transactions"""
        md_db = self._get_db(db_type="md")
        if CURRENT_TRANSACTION in md_db:
            LOG.warning(
                f"Cleaning database {self.get_db_path()} (last transaction failed)"
            )
            data = md_db[CURRENT_TRANSACTION]
            merged_data = data["merged_data"]
            for tax_id in slurm_tqdm(
                merged_data.keys(), desc="cleaning DB", disable=True
            ):
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
        archive = self._archive_stem(file_path)
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
        db_thread.join()
        LOG.debug("Database thread finished")

    def _archive_stem(self, path: str | Path) -> str:
        return Path(path).stem.split(".")[0]

    def _add_data_to_db(self, data):
        with self._db_lock:
            LOG.debug("Locking database transaction")
            md_db = self._get_db(db_type="md")
            md_db[CURRENT_TRANSACTION] = data
            archives = md_db.get(ARCHIVES, [])
            archive = self._archive_stem(data["archive"])

            last_valid_ids = {}

            merged_data = data["merged_data"]
            # warning, tax_id is a str
            LOG.debug(f"Add data to: {self.get_db_path()}")
            for tax_id, tdata in slurm_tqdm(
                merged_data.items(),
                desc=f"Pushing {archive}",
                position=1,
                leave=False,
                disable=True,
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

                for i, (source, counter) in slurm_tqdm(
                    enumerate(zip(sources, counters)),
                    desc="Adding counters and sources",
                    total=len(counters),
                    position=2,
                    leave=False,
                    disable=True,
                ):
                    current_id = last_valid_id + i + 1
                    batch_counters[current_id] = counter
                    batch_sources[current_id] = source

                LOG.debug(
                    f"Starting counters DB transaction with {len(batch_counters)} counters"
                )
                with counter_db.transact():
                    for current_id, counter in slurm_tqdm(
                        batch_counters.items(),
                        desc="Counters",
                        position=3,
                        leave=False,
                        disable=True,
                    ):
                        counter_db[current_id] = counter
                LOG.debug("Ending counters transaction")
                LOG.debug(
                    f"Starting sources DB transaction with {len(batch_sources)} sources"
                )
                with source_db.transact():
                    for current_id, source in slurm_tqdm(
                        batch_sources.items(),
                        desc="Sources",
                        position=3,
                        leave=False,
                        disable=True,
                    ):
                        source_db[current_id] = source
                LOG.debug("Ending sources transaction")

            LOG.debug(f"Data added: {self.get_db_path()}, ending transaction")
            # end transaction
            del md_db[CURRENT_TRANSACTION]
            for tax_id, last_valid_id in last_valid_ids.items():
                self._set_last_valid_id(tax_id=tax_id, last_valid_id=last_valid_id)
            archives.append(archive)
            md_db[ARCHIVES] = archives
            LOG.debug("Transaction ended successfully")
        LOG.debug("Releasing database transaction")

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

    @staticmethod
    def list_dbs(base_path: str | Path) -> list:
        """Get list of available DB (different configs)"""
        base = Path(base_path)
        return [
            [parent.name, child.name]
            for parent in base.iterdir()
            if parent.is_dir()
            for child in parent.iterdir()
            if child.is_dir()
        ]
