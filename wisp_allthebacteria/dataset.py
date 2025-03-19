from collections import defaultdict
from itertools import cycle, product
import logging
import random
from typing import Generator
import numpy as np
import xgboost as xgb
from pathlib import Path
from database import Database
from api import API

LOG = logging.getLogger(__name__)


class Dataset:
    def __init__(self, database: Database, api: API):
        LOG.debug(f"Database({locals()})")
        self._db = database
        self._api = api

        self._column_names = ["".join(p) for p in product("ATGC", repeat=4)]
        self._column_index = {
            name: index for index, name in enumerate(self._column_names)
        }

    def analyse(self, tax_ids: int | list[int]) -> dict:
        if isinstance(tax_ids, int):
            tax_ids = [tax_ids]

        analysis_result = {
            "samples": {},
            "total_samples": 0,
            "rank_counts": defaultdict(
                lambda: defaultdict(lambda: {"count": 0, "tax_id": None})
            ),
        }

        for tax_id in tax_ids:
            # get sample count and additional information
            sample_count = self._db.count(tax_id)
            tax_info = self._api[tax_id]
            scientific_name = tax_info.get("ScientificName", "Unknown")
            rank = tax_info.get("Rank", "Unknown")
            division = tax_info.get("Division", "Unknown")

            # store sample information
            analysis_result["sample"][tax_id] = {
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

        return analysis_result

    def by_rank_generator(
        self, rank: str, batch_size: int, normalize: str | None = None, seed: int = 2025
    ) -> Generator[xgb.DMatrix, None, None]:
        LOG.debug(f"By rank generator for {rank=}, {batch_size=}, {normalize=}")
        tax_ids_by_rank = self._get_tax_ids_by_rank(rank)

        # total samples and weights for each rank_tax_id
        rank_tax_ids = list(tax_ids_by_rank.keys())
        counts = [self._db.count(tax_ids) for tax_ids in tax_ids_by_rank.values()]
        total_samples = sum(counts)
        weights = [count / total_samples for count in counts]

        # generator for each rank_tax_id
        generators = {
            rank_tax_id: self._db.sample_generator(tax_ids, seed=seed, target="counter")
            for rank_tax_id, tax_ids in tax_ids_by_rank.items()
        }

        data_batch = []
        labels_batch = []
        batch_count = 0

        random_instance = random.Random(seed)

        LOG.debug("By rank generator ready - first batch might be slower")
        while generators:
            # select a rank_tax_id based on the calculated weights
            rank_tax_id = random_instance.choices(rank_tax_ids, weights=weights, k=1)[0]
            generator = generators[rank_tax_id]

            try:
                counter = next(generator)
            except StopIteration:
                # remove the generator if it's exhausted
                LOG.debug(f"Generator for {rank_tax_id=} exhausted")
                del generators[rank_tax_id]
                rank_tax_ids.remove(rank_tax_id)
                counts.pop(rank_tax_ids.index(rank_tax_id))
                total_samples = sum(counts)
                weights = [count / total_samples for count in counts]
                continue

            # prepare data and labels for DMatrix
            row = self._counter_to_row(counter=counter, normalize=normalize)
            data_batch.append(row)
            labels_batch.append(rank_tax_id)

            # yield a DMatrix if batch is filled
            if len(data_batch) >= batch_size:
                dmatrix = self._data2DMatrix(data_batch, labels_batch)
                batch_count += 1
                LOG.debug(f"Sending batch {batch_count} - {len(labels_batch)} rows")
                yield dmatrix
                data_batch = []
                labels_batch = []

        # Yield any remaining data as a final DMatrix
        if data_batch:
            dmatrix = self._data2DMatrix(data_batch, labels_batch)
            LOG.debug(f"Sending (last) batch {batch_count} - {len(labels_batch)} rows")
            yield dmatrix

    @staticmethod
    def serialize_dmatrix(dmatrix: xgb.DMatrix, file_path: str | Path) -> None:
        file_path = Path(file_path)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        dmatrix.save_binary(file_path)

    @staticmethod
    def deserialize_dmatrix(file_path: str | Path) -> xgb.DMatrix:
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"{file_path} DMatrix not found.")
        return xgb.DMatrix(file_path)

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

    def _get_tax_ids_by_rank(self, rank: str) -> dict[int, list[int]]:
        rank_mapping = defaultdict(list)
        for tax_id in self._db.get_tax_ids():
            lineage_ex = self._api[tax_id].get("LineageEx", [])
            for entry in lineage_ex:
                if entry["Rank"] == rank:
                    rank_mapping[int(entry["TaxId"])].append(tax_id)
                    break

        return dict(rank_mapping)

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
