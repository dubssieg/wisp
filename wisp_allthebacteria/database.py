import logging
from pathlib import Path
import sys
import threading
import queue
from typing import Literal

from diskcache import Cache
from functools import lru_cache

from reader import Reader
from utils import space_format

LOG = logging.getLogger(__name__)

DB_TYPE = Literal["counter", "source", "index", "md"]
LAST_VALID_ID = "_last_valid_id_"
CURRENT_TRANSACTION = "_current_transaction_"
ARCHIVES = "_archives_"
COMMON = "_common_"


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
        LOG.debug(f"Database({locals()})")
        self._dbs_path = Path(dbs_path).resolve()
        self._kmer_size = kmer_size
        self._reader = reader
        self._window_size = window_size
        self._step = step
        self._full = full
        self.clean()
        self._db_lock = threading.Lock()
        self._task_queue = queue.Queue()
        self._worker_thread = threading.Thread(target=self._db_worker, daemon=True)
        self._worker_thread.start()

    def push_file(self, file_path: str | Path):
        """Add content."""
        file_path = Path(file_path).resolve()
        file_name = file_path.name
        LOG.info(f"pushing file: {file_path}")
        md_db = self._get_db(db_type="md")
        archive = self._archive_stem(file_path)
        archives = md_db.get(ARCHIVES, [])

        if archive in archives:
            LOG.info(f"[{archive}] Already in DB: skip")
            return

        data = self._reader.process_file(
            file_path=file_path,
            kmer_size=self._kmer_size,
            window_size=self._window_size,
            step=self._step,
            full=self._full,
        )
        LOG.info(f"[{file_name}] File processed: data queuing for DB insertion")
        self._task_queue.put(data)
        LOG.debug(f"[{file_name}] Data queued for DB insertion")

    def _db_worker(self):
        """Thread worker managing queue"""
        while True:
            data = self._task_queue.get()  # block if empty
            if data is None:
                break
            try:
                self._add_data_to_db(data)
            except Exception:
                LOG.exception("Database worker error (add_data_to_db)")
                raise
            finally:
                self._task_queue.task_done()

    def _add_data_to_db(self, data):
        LOG.debug("Waiting lock for DB insertion")
        with self._db_lock:
            LOG.debug("Acquiring lock for DB insertion")
            LOG.debug(f"Adding couter and sources : {len(data)} tax_ids")
            md_db = self._get_db(db_type="md")
            md_db[CURRENT_TRANSACTION] = data
            archives = md_db.get(ARCHIVES, [])
            archive = self._archive_stem(data["archive"])

            last_valid_ids = {}

            merged_data = data["merged_data"]
            # warning, tax_id is a str
            LOG.debug(f"Adding data to DB: {self.get_db_path()}")
            for tax_id, tdata in merged_data.items():
                LOG.debug(f"Adding couter and sources : {tax_id}")
                tax_id = self._parse_tax_id(tax_id)
                last_valid_id = self._get_last_valid_id(tax_id)
                counters = tdata["counters"]
                sources = tdata["sources"]
                counter_db = self._get_db(db_type="counter", tax_id=tax_id)
                source_db = self._get_db(db_type="source", tax_id=tax_id)

                # populate and add batch of counters/sources
                batch_counters = {}
                batch_sources = {}

                for i, (source, counter) in enumerate(zip(sources, counters)):
                    current_id = last_valid_id + i + 1
                    batch_counters[current_id] = counter
                    batch_sources[current_id] = source

                LOG.debug(
                    f"Starting counters DB transaction with {len(batch_counters)} counters"
                )
                with counter_db.transact():
                    for current_id, counter in batch_counters.items():
                        counter_db[current_id] = counter
                LOG.debug("Ending counters transaction")
                LOG.debug(
                    f"Starting sources DB transaction with {len(batch_sources)} sources"
                )
                with source_db.transact():
                    for current_id, source in batch_sources.items():
                        source_db[current_id] = source
                LOG.debug("Ending sources transaction")

            LOG.debug("Counters and sources added - ending transaction")
            # end transaction
            del md_db[CURRENT_TRANSACTION]
            for tax_id, last_valid_id in last_valid_ids.items():
                self._set_last_valid_id(tax_id=tax_id, last_valid_id=last_valid_id)
            archives.append(archive)
            md_db[ARCHIVES] = archives
            LOG.debug("Transaction ended successfully, releasing DB lock")
        LOG.debug("DB lock released")

    def wait_for_completion(self):
        """Should be added at the end."""
        LOG.debug("Waiting for DB insertions to complete...")
        self._task_queue.join()
        self.stop_worker()
        LOG.debug("All DB insertions completed")

    def stop_worker(self):
        """Stop worker thread."""
        LOG.debug("Stopping DB worker")
        if self._task_queue.empty():
            self._task_queue.put(None)

        self._worker_thread.join()
        LOG.debug("DB worker stopped")

    def get_counter(self, tax_id: int, num: int) -> dict | None:
        counter_db = self._get_db(db_type="counter", tax_id=tax_id)
        return counter_db.get(num, None)

    def get_source(self, tax_id: int, num: int) -> dict | None:
        source_db = self._get_db(db_type="source", tax_id=tax_id)
        return source_db.get(num, None)

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
            counters[tax_id] = self._get_last_valid_id(tax_id) + 1

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
            for tax_id in merged_data.keys():
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

    def _archive_stem(self, path: str | Path) -> str:
        return Path(path).stem.split(".")[0]

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
        # <base_dbs>/md/common
        # <base_dbs>/md/<tax_id>
        # <base_dbs>/md/counters/<tax_id>
        # <base_dbs>/md/sources/<tax_id>

        path = self.get_db_path()
        path /= db_type
        if tax_id:
            path /= str(tax_id)
        else:
            path /= COMMON

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
