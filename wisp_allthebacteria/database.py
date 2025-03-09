from collections import defaultdict
from itertools import product
from pathlib import Path
import pickle
import random
from typing import Any, Generator
import numpy as np
from xgboost import DMatrix
from tqdm.auto import tqdm
from api import API


class Database:
    def __init__(self, path: str | Path, api: API):
        self._path = Path(path)
        self._api = api
        self._column_names = ["".join(p) for p in product("ATGC", repeat=4)]
        self._column_index = {
            name: index for index, name in enumerate(self._column_names)
        }
        self._dmatrix_filters = dict()

    def get_tax_ids(self) -> list:
        """All DB tax_ids"""
        return [
            int(dir.name)
            for dir in self._path.iterdir()
            if dir.is_dir() and dir.name.isdigit()
        ]

    def get_tax_ids_by_rank(self, rank: str) -> dict[int, list[int]]:
        rank_mapping = defaultdict(list)
        for tax_id in self.get_tax_ids():
            lineage_ex = self._api[tax_id].get("LineageEx", [])
            for entry in lineage_ex:
                if entry["Rank"] == rank:
                    rank_mapping[entry["TaxId"]].append(tax_id)
                    break

        return dict(rank_mapping)

    def _get_tax_id_files(self, tax_id: int) -> list:
        """get all files for a tax_id"""
        tax_id_path = self._path / str(tax_id)
        return list(tax_id_path.glob("*.pkl"))

    def _get_file_data(self, file_path: Path) -> dict:
        with open(file_path, "rb") as file:
            return pickle.load(file)

    # TOFIX: might use to much RAM
    def _get_tax_id_data(self, tax_id: int) -> dict:
        """Read file and get data. Only pickle for now."""
        loaded_data = dict()
        for pickle_file in self._get_tax_id_files(tax_id):
            loaded_data[pickle_file.stem] = self._get_file_data(pickle_file)
        return loaded_data

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

    @staticmethod
    def serialize_dmatrix(dmatrix: DMatrix, file_path: str | Path):
        """Dmatrix serialization."""
        file_path = Path(file_path)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        dmatrix.save_binary(file_path)

    @staticmethod
    def deserialize_dmatrix(file_path: str | Path) -> DMatrix:
        """Dmatrix deserialization."""
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"{file_path} DMatrix not found.")
        return DMatrix(file_path)

    def get_tax_id_classes(self, rank: str) -> int:
        return [int(tax_id) for tax_id in self.get_tax_ids_by_rank(rank).keys()]

    # def add_filter(self, key: str, value: Any):
    #     self._dmatrix_filters[key] = value

    # def get_tax_id_from_filter(self):
    #     tax_ids = list()
    #     for key, value in self._dmatrix_filters.items():
    #         if key == "rank":
    #             tax_ids_by_rank = self.get_tax_ids_by_rank(value)

    def _get_files_by_rank(self, rank: str):
        files_by_rank = dict()
        for rank_tax_id, tax_ids in self.get_tax_ids_by_rank(rank).items():
            files_by_rank[rank_tax_id] = list()
            for tax_id in tax_ids:
                files_by_rank[rank_tax_id].extend(self._get_tax_id_files(tax_id))
        return files_by_rank

    def make_dmatrix(
        self,
        rank: str,
        batch_size: int,
        sample_limit_by_tax_id: int | None = None,
        normalize: str | None = None,
        max_samples: int | None = None,
    ) -> DMatrix | Generator[DMatrix, None, None]:
        """rank can be "phylum, kingdom", etc.
        if batch_size is set, it will create a generator instead of a DMatrix"""
        files_by_rank = self._get_files_by_rank(rank)
        for key in files_by_rank:
            # FIXME: seed
            random.shuffle(files_by_rank[key])

        # by_rank_iterator = {f: 0 for f in files_by_rank.keys()}
        # by_rank_loop_count = {f: 0 for f in files_by_rank.keys()}
        current_batch_count = 0
        total_count = 0

        # get one file for each rank and when all files for a rank have been read, shuffle them
        data = list()
        labels = list()
        can_add_samples = True
        pbar = tqdm(total=max_samples, desc="Building Samples")
        while can_add_samples:
            # clean
            to_remove = list()
            for rank_tax_id, files in files_by_rank.items():
                if not files:
                    to_remove.append(rank_tax_id)
            files_by_rank = {
                k: v for k, v in files_by_rank.items() if k not in to_remove
            }

            # for each rank, take one file read it and remove it from the list
            for rank_tax_id, files in files_by_rank.items():
                selected_file = files.pop()
                content = self._get_file_data(selected_file)
                counters = content["counters"]
                for counter in counters:
                    row = np.zeros(len(self._column_names))
                    if normalize == "sum":
                        counter = self._normalize_sum(counter)
                    elif normalize == "min_max":
                        counter = self._normalize_min_max(counter)
                    # for each "ATGC", etc. count
                    # FIXME: optimize this
                    for col_name, value in counter.items():
                        row[self._column_index[col_name]] = value
                    data.append(row)
                    labels.append(rank_tax_id)
                    total_count += 1
                    pbar.update(1)
                    if max_samples is not None and total_count >= max_samples:
                        can_add_samples = False
                        break
                    current_batch_count += 1
                    if current_batch_count >= batch_size:
                        yield self._data2DMatrix(data, labels, shuffle=True)
                        data = list()
                        labels = list()
                        current_batch_count = 0
        if data:
            yield self._data2DMatrix(data, labels, shuffle=True)
            data = list()
            labels = list()
            current_batch_count = 0

    def _data2DMatrix(self, data: list, labels: list, shuffle: bool = False) -> DMatrix:
        data = np.array(data)
        labels = np.array([int(label) for label in labels])
        if shuffle:
            permuted_id = np.random.permutation(len(labels))
            data = data[permuted_id]
            labels = labels[permuted_id]
        dmatrix = DMatrix(data, label=labels)
        return dmatrix


# TMP
# FIXME: Ugly code but it will do before next refactor
class DMatrixGeneratorFactory:
    def __init__(
        self,
        database: Database,
        rank,
        sample_limit_by_tax_id,
        normalize,
        batch_size,
        max_samples,
    ):
        self.database = database
        self.rank = rank
        self.sample_limit_by_tax_id = sample_limit_by_tax_id
        self.normalize = normalize
        self.batch_size = batch_size
        self.max_samples = max_samples

    def make_dmatrix_generator(self):
        return self.database.make_dmatrix(
            rank=self.rank,
            sample_limit_by_tax_id=self.sample_limit_by_tax_id,
            normalize=self.normalize,
            batch_size=self.batch_size,
            max_samples=self.max_samples,
        )
