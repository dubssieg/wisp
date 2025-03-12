from collections import defaultdict
import hashlib
from itertools import product
from pathlib import Path
import pickle
import random
from typing import Any, Generator
import numpy as np
import xgboost as xgb
from tqdm.auto import tqdm
from api import API

INDEX_VERSION = 1


class Database:
    def __init__(self, path: str | Path, api: API):
        self._path = Path(path)
        self._api = api
        self._column_names = ["".join(p) for p in product("ATGC", repeat=4)]
        self._column_index = {
            name: index for index, name in enumerate(self._column_names)
        }
        self._dmatrix_filters = dict()

    def _get_indexed_data_by_rank(self, rank: str):
        # for each tax_id, read files and index sequences
        index_file = (
            self._path / f"index_by_{rank}_{self._db_signature()}_{INDEX_VERSION}.idx"
        )
        # found: load
        if index_file.exists():
            with open(index_file, "rb") as file:
                return pickle.load(file)
        # not found: compute and save it
        indexed_data = defaultdict(list)
        for rank_tax_id, tax_ids in tqdm(
            self.get_tax_ids_by_rank(rank).items(), desc="Indexing rank", position=1
        ):
            for tax_id in tqdm(
                tax_ids, desc="indexing tax_id", position=2, leave=False
            ):
                tax_id_data = self._get_tax_id_data(tax_id)
                for file_name, data in tqdm(
                    tax_id_data.items(), desc="indexing file", position=3, leave=False
                ):
                    for pos, source in enumerate(data["sources"]):
                        info = self._parse_source(source)
                        info["file"] = file_name
                        info["tax_id"] = tax_id
                        info["pos"] = pos
                        indexed_data[rank_tax_id].append(info)

        with open(index_file, "wb") as file:
            pickle.dump(indexed_data, file)

        return indexed_data

    @staticmethod
    def _parse_source(source: str) -> dict:
        description = source["description"].split()
        return {
            "id": description[0].split(".")[0],
            "len": int(description[1].split("len=")[1]),
        }

    def _db_signature(self):
        """list all files to compute a signature"""
        all_files = set()
        for tax_id in self.get_tax_ids():
            for file in self._get_tax_id_files(tax_id):
                all_files.add(file.stem)
        all_files_str = ",".join(sorted(all_files))
        hash_object = hashlib.sha256()
        hash_object.update(all_files_str.encode())
        return hash_object.hexdigest()

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

    def _get_file_data_by_name_and_tax_id(self, file_name: str, tax_id: int) -> dict:
        file_path = (self._path / str(tax_id) / file_name).with_suffix(".pkl")
        return self._get_file_data(file_path)

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
    def serialize_dmatrix(dmatrix: xgb.DMatrix, file_path: str | Path):
        """Dmatrix serialization."""
        file_path = Path(file_path)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        dmatrix.save_binary(file_path)

    @staticmethod
    def deserialize_dmatrix(file_path: str | Path) -> xgb.DMatrix:
        """Dmatrix deserialization."""
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"{file_path} DMatrix not found.")
        return xgb.DMatrix(file_path)

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

    def make_dmatrix_generator(
        self,
        rank: str,
        batch_size: int,
        sample_limit_by_tax_id: int | None = None,
        normalize: str | None = None,
        max_samples: int | None = None,
    ) -> xgb.DMatrix | Generator[xgb.DMatrix, None, None]:
        """rank can be "phylum, kingdom", etc.
        if batch_size is set, it will create a generator instead of a DMatrix"""
        indexed_data = self._get_indexed_data_by_rank(rank)
        # flatten
        # TODO: balancing strategies
        indexed_data_flat = [
            (rank_tax_id, indexed_item)
            for rank_tax_id, indexed_items in indexed_data.items()
            for indexed_item in indexed_items
        ]
        if max_samples is None or max_samples > len(indexed_data_flat):
            max_samples = len(indexed_data_flat)
        # prepare dataset
        random.shuffle(indexed_data_flat)
        sample_count = 0
        pbar = tqdm(total=max_samples, desc="Making samples", position=1)
        while sample_count < max_samples:
            # extract sample
            batch = [
                indexed_data_flat.pop() for _ in range(batch_size) if indexed_data_flat
            ]
            # sort item by file
            sorted_by_file = defaultdict(list)
            for rank_tax_id, indexed_item in batch:
                file_name = indexed_item["file"]
                tax_id = indexed_item["tax_id"]
                sorted_by_file[(tax_id, file_name)].append((rank_tax_id, indexed_item))
            # read files and fill batch
            data = list()
            labels = list()
            for (tax_id, file_name), sorted_item in tqdm(
                sorted_by_file.items(), desc="Read files", position=2, leave=False
            ):
                file_data = self._get_file_data_by_name_and_tax_id(
                    tax_id=tax_id, file_name=file_name
                )
                for rank_tax_id, indexed_item in tqdm(
                    sorted_item, desc="extract content", position=3, leave=False
                ):
                    pos = indexed_item["pos"]
                    counter = file_data["counters"][pos]
                    data.append(self._counter_to_row(counter, normalize))
                    labels.append(rank_tax_id)
                    sample_count += 1
                    pbar.update()
            yield self._data2DMatrix(data, labels)

    def _counter_to_row(self, counter: dict, normalize: str) -> np.array:
        row = np.zeros(len(self._column_names))
        if normalize == "sum":
            counter = self._normalize_sum(counter)
        elif normalize == "min_max":
            counter = self._normalize_min_max(counter)
        # for each "ATGC", etc. count
        # FIXME: optimize this
        for col_name, value in counter.items():
            row[self._column_index[col_name]] = value
        return row

    def make_dmatrix(
        self,
        rank: str,
        sample_limit_by_tax_id: int | None = None,
        normalize: str | None = None,
    ) -> xgb.DMatrix:
        """rank can be "phylum, kingdom", etc."""
        # sort tax_id in DB by given rank
        tax_ids_by_rank = self.get_tax_ids_by_rank(rank)

        data = []
        labels = []

        # for each rank
        for rank_tax_id, in_rank_tax_ids in tqdm(
            tax_ids_by_rank.items(), position=0, desc="ranks"
        ):
            # for each tax_id in this rank
            for in_rank_tax_id in tqdm(
                in_rank_tax_ids, position=1, leave=False, desc=f"rank {rank_tax_id}"
            ):
                tax_id_data = self._get_tax_id_data(in_rank_tax_id)
                # for each data found in archive (assembly/*tar.gz)
                for file_name, tax_id_data_i in tqdm(
                    tax_id_data.items(),
                    position=2,
                    leave=False,
                    desc=f"sample {in_rank_tax_id}",
                ):
                    sample_count = 0
                    # for each sequence count from *.fa file
                    for counter in tqdm(
                        tax_id_data_i["counters"],
                        position=3,
                        leave=False,
                        desc=f"count for {file_name}",
                    ):
                        sample_count += 1
                        if (
                            sample_limit_by_tax_id is not None
                            and sample_count > sample_limit_by_tax_id
                        ):
                            break
                        # for each "ATGC", etc. count
                        row = self._counter_to_row(counter=counter, normalize=normalize)
                        data.append(row)
                        labels.append(rank_tax_id)

        # convert
        return self._data2DMatrix(data=data, labels=labels, shuffle=True)

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
