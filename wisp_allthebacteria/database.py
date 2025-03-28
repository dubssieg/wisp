import concurrent
import logging
import random
import shutil
import sys
from functools import lru_cache
from pathlib import Path
from typing import Any, Generator, Literal
from tqdm.auto import tqdm

from diskcache import FanoutCache
from fastakmer import FastaKmer
from utils import decompress, space_format, get_weights, sample_count_estimation

LOG = logging.getLogger(__name__)

DB_TYPE = Literal["counter", "index", "md"]
LAST_VALID_IDS = "_last_valid_ids_"
ARCHIVE_TRANSACTION = "_archive_transaction_"
ARCHIVES = "_archives_"


IDX_DB_INFO = "_idx_db_info_"


class Database:
    def __init__(
        self,
        kmer_sizes: int | list[int],
        window_size: int,
        step: int,
        full: bool,
        dbs_path: str | Path,
        fanout_shards: int,
        compression: str | None | bool,
    ):
        LOG.debug(f"Database({locals()})")
        self._dbs_path = Path(dbs_path).resolve()
        self._window_size = window_size
        self._step = step
        self._full = full
        self._fanout_shards = fanout_shards
        self._compression = compression if compression else None
        if isinstance(kmer_sizes, int):
            kmer_sizes = [kmer_sizes]
        self._kmer_sizes = sorted(kmer_sizes)
        self._last_valid_ids = None  # cache

    def get_data(self, tax_id: int, num: int) -> dict | None:
        db = self._get_db(tax_id)
        data = db.get(num, None)
        if data and self._compression:
            return decompress(data, format=self._compression)
        return data

    def count(self, tax_ids: int | str | list[int | str]) -> int:
        """sample count"""
        if isinstance(tax_ids, (int, str)):
            tax_ids = [tax_ids]
        return sum(self._get_last_valid_id(tax_id) + 1 for tax_id in tax_ids)

    def counts(self) -> dict:
        """All counts"""
        return {tid: lvi + 1 for tid, lvi in self._get_last_valid_ids().items()}

    def tax_ids_to_sample_ids(self, tax_ids: int | list[int]) -> list[tuple[int, int]]:
        """Given a tax_id or a list of tax_id, get available samples_ids: DB(tax_id, num)
        Usage: self.get_data(*sample_id)"""
        if isinstance(tax_ids, int):
            tax_ids = [tax_ids]

        sample_ids = []
        for tax_id in tax_ids:
            sample_ids.extend((tax_id, i) for i in range(self.count(tax_id)))
        return sample_ids

    def get_sample_count_estimation(
        self, tax_ids: int | list[int], balance_factor: float = 0.0
    ) -> int:
        if isinstance(tax_ids, int):
            tax_ids = [tax_ids]
        counts = [self.count(tax_id) for tax_id in tax_ids]
        return sample_count_estimation(counts, balance_factor)

    def sample_generator(
        self,
        tax_ids: int | list[int],
        seed: int = 2025,
        target: Literal["counter", "md", "sample_id"] = "counter",
        balance_factor: float = 0.0,
    ) -> Generator[tuple[int, int], None, None]:
        """Generator version of tax_ids_to_sample_ids, with proportional sampling."""
        if isinstance(tax_ids, int):
            tax_ids = [tax_ids]
        random_instance = random.Random(seed)

        LOG.debug(f"Initializing sample generator using {len(tax_ids)} tax_ids")
        counts = [self.count(tax_id) for tax_id in tax_ids]
        weights = get_weights(counts, balance_factor)

        # initialize a list of iterators for each tax_id with random order
        iterators = [
            iter(
                random_instance.sample(
                    range(count),
                    k=count,
                )
            )
            for count in counts
        ]
        LOG.debug("Sample generator initialized")

        while iterators:
            # Select a tax_id based on the calculated weights
            tax_id = random_instance.choices(tax_ids, weights=weights, k=1)[0]
            index = tax_ids.index(tax_id)
            iterator = iterators[index]

            try:
                num = next(iterator)
                sample_id = (tax_id, num)
                match target:
                    case "counter" | "md":
                        yield self.get_data(*sample_id)[target]
                    case "sample_id":
                        yield sample_id
                    case _:
                        raise ValueError(f"{target=}")
            except StopIteration:
                LOG.debug(f"Generator for {tax_id=} exhausted")
                iterators.pop(index)
                tax_ids.pop(index)
                counts.pop(index)
                # Recalculate weights after removing an exhausted tax_id
                weights = get_weights(counts, balance_factor)
        LOG.debug("Sample generator exhausted")

    def get_tax_ids(self) -> list:
        """Get available tax_ids"""
        counters_dir = self.get_db_path() / "counter"
        return [
            int(dir.name)
            for dir in counters_dir.iterdir()
            if dir.is_dir() and dir.name.isdigit()
        ]

    def get_info(self, as_str=False) -> dict:
        """Get DB infos: tax_ids, number of counters, etc."""
        info = self.get_index(IDX_DB_INFO)
        if not info:
            LOG.debug("Computing DB infos...")
            # tax_ids
            tax_ids = self.get_tax_ids()

            # archives pushed
            archives = self.get_archives()

            # counters for each tax_id
            counters = {}
            for tax_id in tqdm(tax_ids, desc="Counting samples..."):
                counters[tax_id] = self.count(tax_id)
            total_samples = sum(counters.values())

            # path
            info = {
                "tax_ids": tax_ids,
                "archives": archives,
                "counters": counters,
                "total_samples": total_samples,
            }
            LOG.debug("DB infos - DONE")
            self.set_index(IDX_DB_INFO, info)

        if as_str:
            info_lines = []
            info_lines.append("\n=== Paths/Config ===")
            info_lines.append(f"Current base path: {self.get_db_path()}")
            info_lines.append("\n=== Sequences ===")
            info_lines.append(
                f"{len(info['tax_ids'])} tax_ids -> [tax_id] sample_count:"
            )
            items_per_line = 5
            counters_str = {
                tax_id: space_format(num)
                for tax_id, num in sorted(info["counters"].items())
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
            info_lines.append(f"Total: {space_format(info['total_samples'])} samples")

            info_lines.append(f"\n==== Archives pushed ({len(info['archives'])}) ====")
            info_lines.append(", ".join(sorted(info["archives"])))

            return "\n".join(info_lines)

        return info

    def set_index(self, key: Any, value: Any) -> None:
        """Everything indexed will be deleted next time a file is pushed on DB"""
        LOG.debug(f"Add/update INDEX DB for {key=}")
        in_db = self._get_db(db_type="index")
        in_db[key] = value

    def get_index(self, key: Any) -> Any:
        """Indexed value or None"""
        in_db = self._get_db(db_type="index")
        return in_db.get(key, None)

    def clear_index(self) -> None:
        LOG.debug("Clearing INDEX DB...")
        idx_db = self._get_db(db_type="index")
        idx_db.clear()
        LOG.debug("INDEX DB is now empty")

    def get_archives(self) -> list[str]:
        """Pushed archives"""
        md_db = self._get_db(db_type="md")
        return md_db.get(ARCHIVES, [])

    def has_archive(self, path: str | Path) -> bool:
        """Check if archive already was pushed in base."""
        return self._archive_stem(path) in self.get_archives()

    def _archive_stem(self, path: str | Path) -> str:
        return Path(path).stem.split(".")[0]

    def _get_last_valid_ids(self) -> dict:
        """Get last inserted ids"""
        if self._last_valid_ids is None:
            md_db = self._get_db(db_type="md")
            self._last_valid_ids = md_db.get(LAST_VALID_IDS, {})
        return self._last_valid_ids

    def _get_last_valid_id(self, tax_id: int) -> int:
        return self._get_last_valid_ids().get(tax_id, -1)

    def _parse_tax_id(self, tax_id: str | int) -> str | int:
        try:
            return int(tax_id)
        except (ValueError, TypeError):
            return tax_id

    @lru_cache(maxsize=1)
    def _get_db(
        self, tax_id: int | str | None = None, db_type: DB_TYPE = "counter"
    ) -> FanoutCache:
        """Get sub DB"""
        path = self.get_db_path()
        path /= db_type
        if db_type == "counter":
            if tax_id:
                path /= str(tax_id)
            else:
                raise ValueError("tax_id needed")

        path.mkdir(parents=True, exist_ok=True)
        return FanoutCache(path, size_limit=sys.maxsize, shards=self._fanout_shards)

    def get_db_path(self):
        dir_name = f"kmr_{'_'.join(map(str, self._kmer_sizes))}"
        if self._full:
            dir_name += "__full"
        else:
            dir_name += f"__wsz_{self._window_size}__stp_{self._step}"

        if self._compression:
            dir_name += f"__comp_{self._compression}"

        dir_name += f"__shd_{self._fanout_shards}"

        return self._dbs_path / dir_name


class DatabaseBuilder(Database):
    def __init__(
        self,
        kmer_sizes: int | list[int],
        window_size: int,
        step: int,
        full: bool,
        dbs_path: str | Path,
        fanout_shards: int,
        fasta_kmer: FastaKmer,
        insert_threads: int,
        compression: str | bool | None,
        fasta_batch_size: int | None = None,
        merged_data_as_db: bool = False,
        species_count_limit: int | None = None,
    ):
        LOG.debug(f"DatabaseBuilder({locals()})")
        super().__init__(
            kmer_sizes=kmer_sizes,
            window_size=window_size,
            step=step,
            full=full,
            fanout_shards=fanout_shards,
            dbs_path=dbs_path,
            compression=compression,
        )
        self._fasta_kmer = fasta_kmer
        self._fasta_batch_size = fasta_batch_size
        self._insert_threads = insert_threads
        self._merged_data_as_db = merged_data_as_db
        self._species_count_limit = species_count_limit
        self._cancel_transaction()  # if needed

    def push_archive(self, archive_path: str | Path):
        """Add content."""
        archive_path = Path(archive_path).resolve()
        LOG.info(f"Pushing file: {archive_path.name}")

        self._start_transaction()
        for data in self._fasta_kmer.process_archive(
            archive_path=archive_path,
            kmer_sizes=self._kmer_sizes,
            window_size=self._window_size,
            step=self._step,
            full=self._full,
            batch_size=self._fasta_batch_size,
            compression=self._compression,
            merged_data_as_db=self._merged_data_as_db,
            max_count=self._species_count_limit,
            current_counts=self.counts(),
        ):
            if data:
                self._push_merged_data(data)

        self.clear_index()
        self._end_transaction(archive_path)

    def _apply_count_limit(self, data: dict) -> dict:
        """deprecated?"""
        if self._species_count_limit:
            new_merged_data = {}
            for tid, counters in data["merged_data"].items():
                current_count = self.count(tid)
                new_count_limit = self._species_count_limit - current_count
                if new_count_limit > 0:
                    new_merged_data[tid] = counters[:new_count_limit]
                    if len(counters) > new_count_limit:
                        LOG.debug(
                            f"Filtered {len(counters) - new_count_limit} counters for specie {tid}"
                        )
                else:
                    LOG.debug(
                        f"All {len(counters)} counters filtered for specie {tid} due to count limit"
                    )

            data["merged_data"] = new_merged_data

        return data

    def _push_merged_data(self, data) -> dict[int, int]:
        merged_data = data["merged_data"]
        tmp_dir = data["tmp_dir"]
        LOG.debug(f"Adding counters to DB for {len(merged_data)} tax_id")

        last_valid_ids = {}
        insert_counter = 0
        is_db = data["tmp_dir"] is not None

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
                    LOG.exception(f"Pushing counters to DB for tax_id {tax_id}")
                    raise
        if is_db:
            LOG.debug(f"Deleting temporary directory {tmp_dir} ...")
            shutil.rmtree(tmp_dir)
            LOG.debug("Temporary files deleted")

        self._set_last_valid_ids(last_valid_ids)

    def _push_batch(self, tax_id: int, tdata: dict, is_db: bool):
        LOG.debug(f"Processing tax_id: {tax_id}")
        last_valid_id = self._get_last_valid_id(tax_id)
        dst_db = self._get_db(db_type="counter", tax_id=tax_id)

        if is_db:
            src_db = tdata["db"]
            src_last_id = tdata["last_id"]
            LOG.debug(
                f"Starting DB transaction with {src_last_id} counters for tax_id {tax_id}"
            )

            with dst_db.transact():
                for j in range(src_last_id):
                    dst_db[last_valid_id + j + 1] = src_db[j]

            new_last_valid_id = last_valid_id + src_last_id

        else:
            batch_counters = {
                last_valid_id + j + 1: counter for j, counter in enumerate(tdata)
            }

            LOG.debug(
                f"Starting DB transaction : counte with {len(batch_counters)} counters for tax_id {tax_id}"
            )

            with dst_db.transact():
                for current_id, counter in batch_counters.items():
                    dst_db[current_id] = counter

            new_last_valid_id = last_valid_id + len(batch_counters)

        LOG.debug(f"Finished specie: {tax_id}")
        return tax_id, new_last_valid_id

    def _set_last_valid_id(self, tax_id: int | str, last_valid_id: int) -> None:
        """Set last inserted id
        WARNING: not thread safe (self._last_valid_ids)"""
        last_valid_ids = self._get_last_valid_ids()
        last_valid_ids[tax_id] = last_valid_id
        md_db = self._get_db(db_type="md")
        md_db[LAST_VALID_IDS] = last_valid_ids
        self._last_valid_ids = None

    def _set_last_valid_ids(self, updated_last_valid_ids: dict[int | str, int]) -> None:
        """Update last inserted ids
        WARNING: not thread safe (self._last_valid_ids)"""
        last_valid_ids = self._get_last_valid_ids()
        last_valid_ids.update(updated_last_valid_ids)
        md_db = self._get_db(db_type="md")
        md_db[LAST_VALID_IDS] = last_valid_ids
        self._last_valid_ids = None

    def _start_transaction(self) -> None:
        LOG.debug("starting transaction")
        md_db = self._get_db(db_type="md")
        md_db[ARCHIVE_TRANSACTION] = self._get_last_valid_ids()

    def _end_transaction(self, archive: str | Path) -> None:
        LOG.debug(f"ending transaction for {archive}")
        md_db = self._get_db(db_type="md")
        archives: list = md_db.get(ARCHIVES, [])
        archive = self._archive_stem(archive)
        archives.append(archive)
        md_db[ARCHIVES] = archives
        del md_db[ARCHIVE_TRANSACTION]

    def _cancel_transaction(self):
        md_db = self._get_db(db_type="md")
        if ARCHIVE_TRANSACTION in md_db:
            LOG.warning("Cancel transaction - Restoring...")
            last_valid_ids = md_db[ARCHIVE_TRANSACTION]
            self._set_last_valid_ids(last_valid_ids)
            del md_db[ARCHIVE_TRANSACTION]
            LOG.warning("Transaction canceled")
