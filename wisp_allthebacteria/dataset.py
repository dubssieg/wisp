from collections import defaultdict
import concurrent.futures
from itertools import product
import logging
import queue
import random
import threading
import time
from typing import Generator
import concurrent
import numpy as np
import xgboost as xgb
from database import Database
from api import API
from utils import hash, FunctionLogger

LOG = logging.getLogger(__name__)

IDX_ANALYSIS = "_analysis_"
IDX_TID_BY_RANK = "_tid_by_rank_"


class Dataset:
    def __init__(self, database: Database, api: API):
        LOG.debug(f"Database({locals()})")
        self._db = database
        self._api = api

        self._column_names = ["".join(p) for p in product("ATGC", repeat=4)]
        self._column_index = {
            name: index for index, name in enumerate(self._column_names)
        }

    def analyse(self, tax_ids: int | list[int] | None = None) -> dict:
        if tax_ids is None:
            tax_ids = self._db.get_tax_ids()
        if isinstance(tax_ids, int):
            tax_ids = [tax_ids]

        idx_key = (IDX_ANALYSIS, hash(tax_ids))
        analysis_result = self._db.get_index(idx_key)
        if not analysis_result:
            LOG.debug(f"Analysing {len(tax_ids)} tax_ids")
            analysis_result = {
                "samples": {},
                "total_samples": 0,
                "rank_counts": defaultdict(
                    lambda: defaultdict(lambda: {"count": 0, "tax_id": None})
                ),
            }

            for tax_id in tax_ids:
                db_info = self._db.get_info()
                # get sample count and additional information
                sample_count = db_info["counters"][tax_id]
                tax_info = self._api[tax_id]
                scientific_name = tax_info.get("ScientificName", "Unknown")
                rank = tax_info.get("Rank", "Unknown")
                division = tax_info.get("Division", "Unknown")

                # store sample information
                analysis_result["samples"][tax_id] = {
                    "scientific_name": scientific_name,
                    "rank": rank,
                    "division": division,
                    "count": sample_count,
                }
                analysis_result["total_samples"] += sample_count

                # lineage information (ranks)
                lineage_ex = tax_info.get("LineageEx", [])

                # count the number of elements for each parent rank
                for entry in lineage_ex:
                    parent_rank = entry["Rank"]
                    parent_scientific_name = entry["ScientificName"]
                    parent_tax_id = int(entry["TaxId"])

                    # update the rank_counts structure
                    analysis_result["rank_counts"][parent_rank][parent_scientific_name][
                        "count"
                    ] += sample_count
                    analysis_result["rank_counts"][parent_rank][parent_scientific_name][
                        "tax_id"
                    ] = parent_tax_id

            analysis_result["rank_counts"] = {
                k: dict(v) for k, v in analysis_result["rank_counts"].items()
            }
            self._db.set_index(idx_key, analysis_result)

        return analysis_result

    def total_samples(self) -> int:
        db_info = self._db.get_info()
        return db_info["total_samples"]

    def labels(self, rank: str) -> list[str]:
        """All labels in this DB, for this rank"""
        return list(self._get_tax_ids_by_rank(rank).keys())

    def _get_tax_ids_by_rank(self, rank: str) -> dict[int, list[int]]:
        idx_key = (IDX_TID_BY_RANK, rank)
        rank_mapping = self._db.get_index(idx_key)
        if not rank_mapping:
            rank_mapping = defaultdict(list)
            for tax_id in self._db.get_tax_ids():
                lineage_ex = self._api[tax_id].get("LineageEx", [])
                for entry in lineage_ex:
                    if entry["Rank"] == rank:
                        rank_mapping[int(entry["TaxId"])].append(tax_id)
                        break
            rank_mapping = dict(rank_mapping)
            self._db.set_index(idx_key, rank_mapping)

        return rank_mapping


class ByRankGenerator(Dataset):
    def __init__(
        self,
        database: Database,
        api: API,
        rank: str,
        batch_size: int,
        normalize: str | None = None,
        seed: int = 2025,
        buffer_threads: int = 10,
    ):
        super().__init__(database=database, api=api)
        LOG.debug(f"ByRankGenerator for {rank=}, {batch_size=}, {normalize=}")
        self._rank = rank
        self._batch_size = batch_size
        self._batch_count = 0
        self._normalize = normalize
        self._seed = seed
        self._buffer_threads = buffer_threads
        self._lock = threading.Lock()
        self._tax_ids_by_rank = self._get_tax_ids_by_rank(rank)
        self._rank_tids = list(self._tax_ids_by_rank.keys())
        self._max_buffer_size = max(32, batch_size // len(self._rank_tids))
        self._buffers = {rank_tid: queue.Queue() for rank_tid in self._rank_tids}
        self._exhausted = {rank_tid: False for rank_tid in self._rank_tids}
        self._filling = {rank_tid: False for rank_tid in self._rank_tids}
        self._generators = {
            rank_tid: self._db.sample_generator(
                tax_ids, seed=self._seed, target="counter"
            )
            for rank_tid, tax_ids in self._tax_ids_by_rank.items()
        }
        self._data_batch = []
        self._labels_batch = []

        db_info = self._db.get_info()
        self._counts = [
            sum(db_info["counters"][tax_id] for tax_id in tax_ids)
            for tax_ids in self._tax_ids_by_rank.values()
        ]
        self._total_samples = sum(self._counts)
        self._weights = [count / self._total_samples for count in self._counts]
        self._max_batch_count = self.available_batches_count()

        # start filling buffers
        self._executor = concurrent.futures.ThreadPoolExecutor(
            max_workers=self._buffer_threads
        )
        self._terminating = False

    def start(self) -> None:
        LOG.debug("Start filling buffers...")
        self._start_filling_buffers()

    def stop(self) -> None:
        LOG.debug("Stopping threads")
        self._terminating = True
        self._executor.shutdown(wait=True)

    def get(self) -> Generator[xgb.DMatrix, None, None]:

        random_instance = random.Random(self._seed)

        with (
            FunctionLogger(30, self._batch_info),
            FunctionLogger(30, self._buffers_info),
        ):

            while not (
                self._all_generator_done()
                and self._all_buffers_done()
                and self._terminating
            ):
                rank_tid = random_instance.choices(
                    self._rank_tids, weights=self._weights, k=1
                )[0]
                try:
                    counter = self._buffers[rank_tid].get()
                except queue.Empty:
                    if self._is_done(rank_tid):
                        LOG.debug(f"{self._rank} {rank_tid} is DONE - Cleaning")
                        self._remove(rank_tid)
                        continue

                # prepare data and labels for DMatrix
                row = self._counter_to_row(counter=counter, normalize=self._normalize)
                self._data_batch.append(row)
                self._labels_batch.append(rank_tid)

                # yield a DMatrix if batch is filled
                if len(self._data_batch) >= self._batch_size:
                    dmatrix = self._data2DMatrix(self._data_batch, self._labels_batch)
                    self._batch_count += 1
                    LOG.debug(
                        f"Sending batch {self._batch_count} / (max: {self._max_batch_count}) - {len(self._labels_batch)} samples"
                    )
                    self._data_batch = []
                    self._labels_batch = []
                    yield dmatrix

            # yield any remaining data as a final DMatrix
            if self._data_batch:
                self._batch_count += 1
                dmatrix = self._data2DMatrix(self._data_batch, self._labels_batch)
                LOG.debug(
                    f"Sending (last) batch {self._batch_count} - {len(self._labels_batch)} samples"
                )
                yield dmatrix

    def _buffers_info(self) -> str:
        buffer_info = []
        for tid in self._rank_tids:
            qsize = str(self._buffers[tid].qsize())
            filling = self._filling[tid]
            exhausted = self._exhausted[tid]

            if filling:
                qsize = f"^{qsize}^"
            elif exhausted:
                qsize = f"[{qsize}]"

            buffer_info.append(qsize)

        return (
            f"BUFFERS (max: {self._max_buffer_size}, th: {self._buffer_threads}): "
            + "|".join(buffer_info)
        )

    def _batch_info(self) -> None:
        perc = 100 * len(self._labels_batch) / self._batch_size
        return f"BATCH {self._batch_count + 1}: {len(self._labels_batch) / self._batch_size} ({perc:.1f} %)"

    def available_batches_count(self) -> int:
        """How many batches are available"""
        total_samples = self.total_samples()
        return (total_samples + self._batch_size - 1) // self._batch_size

    def _fill_buffers(self):
        """Thread worker dynamically picking the least filled buffer."""
        while not self._all_generator_done() and not self._terminating:
            rank_tid = self._select_next_buffer()
            if rank_tid is None:
                time.sleep(1)
                continue

            try:
                sample = next(self._generators[rank_tid])
                self._buffers[rank_tid].put(sample)
                # self.tmp_log()  # TODO: remove
            except StopIteration:
                self._mark_exhausted(rank_tid)
            self._stop_filling(rank_tid)

    def _start_filling_buffers(self):
        for _ in range(self._executor._max_workers):
            self._executor.submit(self._fill_buffers)

    def _select_next_buffer(self) -> int | None:
        with self._lock:
            if rank_tid := min(
                (
                    rank_tid
                    for rank_tid, q in self._buffers.items()
                    if not self._filling[rank_tid] and q.qsize() < self._max_buffer_size
                ),
                key=lambda rank_tid: self._buffers[rank_tid].qsize(),
                default=None,
            ):
                self._filling[rank_tid] = True
                return rank_tid

    def _all_generator_done(self) -> bool:
        with self._lock:
            return all(self._exhausted[rank_tid] for rank_tid in self._rank_tids)

    def _all_buffers_done(self) -> bool:
        with self._lock:
            return all(
                self._buffers[rank_tid].qsize() == 0 for rank_tid in self._rank_tids
            )

    def _is_rank_done(self, rank_tid: int) -> bool:
        with self._lock:
            return (
                self._buffers[rank_tid].qsize() == 0
                and self._exhausted[rank_tid]
                and not self._filling[rank_tid]
            )

    def _mark_exhausted(self, rank_tid: int):
        LOG.debug(f"Generator for {rank_tid} is exhausted")
        with self._lock:
            self._exhausted[rank_tid] = True

    def _stop_filling(self, rank_tid: int):
        with self._lock:
            self._filling[rank_tid] = False

    def _remove(self, rank_tid: int):
        with self._lock:
            del self._generators[rank_tid]
            self._rank_tids.remove(rank_tid)
            self._counts.pop(self._rank_tids.index(rank_tid))
            self._total_samples = sum(self._counts)
            self._weights = [count / self._total_samples for count in self._counts]

    @staticmethod
    def _normalize_sum(count_dict: dict) -> dict:
        total_count = sum(count_dict.values())
        return {key: value / total_count for key, value in count_dict.items()}

    @staticmethod
    def _normalize_min_max(count_dict) -> dict:
        min_value = min(count_dict.values())
        max_value = max(count_dict.values())
        if max_value == min_value:
            return {key: 1.0 for key in count_dict}
        return {
            key: (value - min_value) / (max_value - min_value)
            for key, value in count_dict.items()
        }

    def _counter_to_row(self, counter: dict, normalize: str) -> np.array:
        row = np.zeros(len(self._column_names))
        if normalize == "sum":
            counter = self._normalize_sum(counter)
        elif normalize == "min_max":
            counter = self._normalize_min_max(counter)
        # for each "ATGC", etc. count
        indices = np.array(
            [self._column_index[col_name] for col_name in counter.keys()]
        )
        values = np.array(list(counter.values()))
        np.put(row, indices, values)
        return row

    def _data2DMatrix(
        self, data: list, labels: list, shuffle: bool = False
    ) -> xgb.DMatrix:
        data = np.array(data)
        labels = np.array([int(label) for label in labels])
        if shuffle:
            permuted_id = np.random.permutation(len(labels))
            data = data[permuted_id]
            labels = labels[permuted_id]
        dmatrix = xgb.DMatrix(data, label=labels)
        return dmatrix
