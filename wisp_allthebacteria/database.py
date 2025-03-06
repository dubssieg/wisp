from collections import defaultdict
from pathlib import Path

from api import API


class Database:
    def __init__(self, path: str | Path, api: API):
        self._path = Path(path)
        self._api = api

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
