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

    # def make_dmatrix(self, rank: str) -> xgb.DMatrix:
    #     tax_ids_by_rank = self.get_tax_ids_by_rank(rank)
    #     for rank_tax_id, samples_tax_id in tax_ids_by_rank.items():
    #         tax_id_data = self._get_tax_id_data(samples_tax_id)

    def make_dmatrix(
        self, rank: str, sample_limit_by_tax_id: int | None = None
    ) -> xgb.DMatrix:
        tax_ids_by_rank = self.get_tax_ids_by_rank(rank)

        data = []
        labels = []

        for rank_tax_id, in_rank_tax_ids in tqdm(
            tax_ids_by_rank.items(), position=0, desc="ranks"
        ):

            for in_rank_tax_id in tqdm(
                in_rank_tax_ids, position=1, leave=False, desc=f"rank {rank_tax_id}"
            ):

                tax_id_data = self._get_tax_id_data(in_rank_tax_id)

                for file_name, tax_id_data_i in tqdm(
                    tax_id_data.items(),
                    position=2,
                    leave=False,
                    desc=f"sample {in_rank_tax_id}",
                ):

                    sample_count = 0
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
                        for col_name, value in counter.items():
                            row[self._column_index[col_name]] = value
                        data.append(row)
                        labels.append(rank_tax_id)

        # Convert data and labels to numpy arrays
        data = np.array(data)
        labels = np.array([int(label) for label in labels])

        # Create DMatrix
        dmatrix = xgb.DMatrix(data, label=labels)
        return dmatrix
