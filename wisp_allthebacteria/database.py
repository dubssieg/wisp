import concurrent
import logging
import shutil
import sys
from functools import lru_cache
from pathlib import Path
from typing import Generator, Literal

from diskcache import FanoutCache
from reader import Reader
from utils import decompress, space_format

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
        fanout_shards: int,
        compressed: bool,
    ):
        LOG.debug(f"Database({locals()})")
        self._dbs_path = Path(dbs_path).resolve()
        self._kmer_size = kmer_size
        self._window_size = window_size
        self._step = step
        self._full = full
        self._fanout_shards = fanout_shards
        self._compressed = compressed

    def get_data(
        self, tax_id: int, num: int, target: Literal["counter", "source"] = "counter"
    ) -> dict | None:
        db = self._get_db(db_type=target, tax_id=tax_id)
        data = db.get(num, None)
        if data and self._compressed:
            return decompress(data)
        return data

    def tax_ids_to_sample_ids(self, tax_ids: int | list[int]) -> list[tuple[int, int]]:
        """Given a tax_id or a list of tax_id, get available samples_ids: DB(tax_id, num)
        Usage: self.get_data(*sample_id)"""
        if isinstance(tax_ids, int):
            tax_ids = [tax_ids]

        sample_ids = []
        for tax_id in tax_ids:
            self._get_last_valid_id(tax_id)
            sample_ids.extend(
                (tax_id, i) for i in range(self._get_last_valid_id(tax_id) + 1)
            )
        return sample_ids

    def tax_ids_to_sample_ids_generator(
        self, tax_ids: int | list[int]
    ) -> Generator[tuple[int, int], None, None]:
        """Generator version of tax_ids_to_sample_ids, looping on tax_ids.
        Take one of tax_id 1, take one, of tax_id 2, etc."""
        if isinstance(tax_ids, int):
            tax_ids = [tax_ids]

        # initialize a list of iterators for each tax_id
        iterators = [
            iter(range(self._get_last_valid_id(tax_id) + 1)) for tax_id in tax_ids
        ]

        while iterators:
            for tax_id, iterator in zip(tax_ids, iterators):
                try:
                    num = next(iterator)
                    yield (tax_id, num)
                except StopIteration:
                    # remove the iterator if it's exhausted
                    iterators.remove(iterator)
                    tax_ids.remove(tax_id)

    def get_tax_ids(self) -> list:
        """Get available tax_ids - TODO: in MD DB"""
        counters_dir = self.get_db_path() / "counter"
        return [
            int(dir.name)
            for dir in counters_dir.iterdir()
            if dir.is_dir() and dir.name.isdigit()
        ]

    def get_info(self, as_str=False) -> dict:
        """Get DB infos: tax_ids, number of counters, etc."""
        # tax_ids
        tax_ids = self.get_tax_ids()

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
            info_lines.append("\n=== Sequences ===")
            info_lines.append(f"{len(tax_ids)} tax_ids -> [tax_id] sample_count:")
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

    def _archive_stem(self, path: str | Path) -> str:
        return Path(path).stem.split(".")[0]

    def _get_last_valid_id(self, tax_id: int) -> int:
        """Get last inserted id"""
        md_db = self._get_db(db_type="md", tax_id=tax_id)
        return md_db.get(LAST_VALID_ID, -1)

    def _parse_tax_id(self, tax_id: str | int) -> str | int:
        try:
            return int(tax_id)
        except (ValueError, TypeError):
            return tax_id

    @lru_cache(maxsize=1)
    def _get_db(self, db_type: DB_TYPE, tax_id: int | None = None) -> FanoutCache:
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

        return FanoutCache(path, size_limit=sys.maxsize, shards=self._fanout_shards)

    def get_db_path(self):

        dir_name = f"km_{self._kmer_size}"
        if self._full:
            dir_name += "__full"
        else:
            dir_name += f"__ws_{self._window_size}__st_{self._step}"

        if self._compressed:
            dir_name += "__comp"

        dir_name += f"__sh_{self._fanout_shards}"

        return self._dbs_path / dir_name


class DatabaseBuilder(Database):
    def __init__(
        self,
        kmer_size: int,
        window_size: int,
        step: int,
        full: bool,
        dbs_path: str | Path,
        fanout_shards: int,
        reader: Reader,
        insert_threads: int,
        compressed: bool,
        fasta_batch_size: int | None = None,
        merged_data_as_db: bool = False,
    ):
        LOG.debug(f"DatabaseBuilder({locals()})")
        super().__init__(
            kmer_size=kmer_size,
            window_size=window_size,
            step=step,
            full=full,
            fanout_shards=fanout_shards,
            dbs_path=dbs_path,
            compressed=compressed,
        )
        self._reader = reader
        self._fasta_batch_size = fasta_batch_size
        self._insert_threads = insert_threads
        self._merged_data_as_db = merged_data_as_db
        self.clean()

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

    def push_file(self, file_path: str | Path):
        """Add content."""
        file_path = Path(file_path).resolve()
        LOG.info(f"Pushing file: {file_path.name}")

        data = self._reader.process_file(
            file_path=file_path,
            kmer_size=self._kmer_size,
            window_size=self._window_size,
            step=self._step,
            full=self._full,
            batch_size=self._fasta_batch_size,
            compressed=self._compressed,
            merged_data_as_db=self._merged_data_as_db,
        )

        self._push_merged_data(data)

    def _push_merged_data(self, data):
        LOG.debug("Adding counter & sources - get DB metadata")
        md_db = self._get_db(db_type="md")
        md_db[CURRENT_TRANSACTION] = data
        archives = md_db.get(ARCHIVES, [])
        archive = self._archive_stem(data["archive"])

        # warning, tax_id is a str
        merged_data = data["merged_data"]
        LOG.debug(
            f"Adding counters & sources to DB for {len(merged_data)} tax_id(s): {self.get_db_path()}"
        )

        last_valid_ids = {}
        insert_counter = 0
        is_db = data["tmp_dir"] is not None

        # 1 - update counters & sources
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self._insert_threads
        ) as executor:
            future_to_taxid = {
                executor.submit(
                    self._push_batch, tax_id=tax_id, tdata=tdata, is_db=is_db
                ): tax_id
                for tax_id, tdata in merged_data.items()
            }

            for future in concurrent.futures.as_completed(future_to_taxid):
                try:
                    tax_id, last_valid_id = future.result()
                    last_valid_ids[tax_id] = last_valid_id
                    insert_counter += 1
                    LOG.debug(f"DB insertions: {insert_counter} / {len(merged_data)}")
                except Exception:
                    LOG.exception(
                        f"Pushing counters & sources to DB for tax_id {tax_id}"
                    )
                    raise

        LOG.debug("All counters & sources added - ending transaction")
        # end transaction
        del md_db[CURRENT_TRANSACTION]

        # 2 - update MD
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self._insert_threads
        ) as executor:
            future_to_taxid = {
                executor.submit(
                    self._set_last_valid_id, tax_id=tax_id, last_valid_id=last_valid_id
                ): tax_id
                for tax_id, last_valid_id in last_valid_ids.items()
            }

            for future in concurrent.futures.as_completed(future_to_taxid):
                tax_id = future_to_taxid[future]
                try:
                    future.result()
                except Exception:
                    LOG.exception(f"Error updating last valid ID for tax_id {tax_id}")
                    raise

        archives.append(archive)
        md_db[ARCHIVES] = archives
        LOG.debug("Transaction ended successfully")

        if is_db:
            LOG.debug(f"Deleting temporary directory {data['tmp_dir']} ...")
            shutil.rmtree(data["tmp_dir"])
            LOG.debug("Temporary files deleted")

    def _push_batch(self, tax_id: int, tdata: dict, is_db: bool):
        LOG.debug(f"Processing tax_id: {tax_id}")
        tax_id = self._parse_tax_id(tax_id)
        last_valid_id = self._get_last_valid_id(tax_id)
        counters = tdata["counters"]
        sources = tdata["sources"]
        counter_db = self._get_db(db_type="counter", tax_id=tax_id)
        source_db = self._get_db(db_type="source", tax_id=tax_id)

        if is_db:
            merged_data_last_id = tdata["last_id"]
            LOG.debug(
                f"Starting DB transactions with {merged_data_last_id} counters & sources for tax_id {tax_id}"
            )

            with counter_db.transact():
                for j in range(merged_data_last_id):
                    counter_db[last_valid_id + j + 1] = counters[j]
            with source_db.transact():
                for j in range(merged_data_last_id):
                    source_db[last_valid_id + j + 1] = sources[j]

            new_last_valid_id = last_valid_id + merged_data_last_id

        else:
            batch_counters = {
                last_valid_id + j + 1: counter for j, counter in enumerate(counters)
            }
            batch_sources = {
                last_valid_id + j + 1: source for j, source in enumerate(sources)
            }

            LOG.debug(
                f"Starting DB transactions with {len(batch_counters)} counters & sources for tax_id {tax_id}"
            )

            with counter_db.transact():
                for current_id, counter in batch_counters.items():
                    counter_db[current_id] = counter
            with source_db.transact():
                for current_id, source in batch_sources.items():
                    source_db[current_id] = source

            new_last_valid_id = last_valid_id + len(batch_counters)

        LOG.debug(f"Finished tax_id: {tax_id}")
        return tax_id, new_last_valid_id

    def _set_last_valid_id(self, tax_id: int, last_valid_id: int) -> None:
        """Get last inserted id"""
        md_db = self._get_db(db_type="md", tax_id=tax_id)
        md_db[LAST_VALID_ID] = last_valid_id
