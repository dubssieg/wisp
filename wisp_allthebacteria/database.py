from collections import defaultdict
from itertools import product
from pathlib import Path
import pickle
import numpy as np
import xgboost as xgb
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

    def _get_tax_id_data(self, tax_id: int) -> dict:
        """Read file and get data. Only pickle for now."""
        tax_id_path = self._path / str(tax_id)
        pickle_files = list(tax_id_path.glob("*.pkl"))
        loaded_data = dict()
        for pickle_file in pickle_files:
            with open(pickle_file, "rb") as file:
                loaded_data[pickle_file.stem] = pickle.load(file)
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
                        row = np.zeros(len(self._column_names))
                        sample_count += 1
                        if (
                            sample_limit_by_tax_id is not None
                            and sample_count > sample_limit_by_tax_id
                        ):
                            break
                        # for each "ATGC", etc. count
                        if normalize == "sum":
                            counter = self._normalize_sum(counter)
                        elif normalize == "min_max":
                            counter = self._normalize_min_max(counter)

                        for col_name, value in counter.items():
                            row[self._column_index[col_name]] = value
                        data.append(row)
                        labels.append(rank_tax_id)

        # convert
        data = np.array(data)
        labels = np.array([int(label) for label in labels])
        dmatrix = xgb.DMatrix(data, label=labels)
        return dmatrix
