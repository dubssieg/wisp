from collections import defaultdict
import logging
from database import Database
from api import API

LOG = logging.getLogger(__name__)


class Dataset:
    def __init__(self, database: Database, api: API):
        LOG.debug(f"Database({locals()})")
        self._db = database
        self._api = api

    def get_tax_ids_by_rank(self, rank: str) -> dict[int, list[int]]:
        rank_mapping = defaultdict(list)
        for tax_id in self._db.get_tax_ids():
            lineage_ex = self._api[tax_id].get("LineageEx", [])
            for entry in lineage_ex:
                if entry["Rank"] == rank:
                    rank_mapping[entry["TaxId"]].append(tax_id)
                    break

        return dict(rank_mapping)
