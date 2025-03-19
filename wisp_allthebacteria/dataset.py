from collections import defaultdict
from itertools import cycle, product
import logging
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


    def analyse_ranks(
        self, tax_ids: int | list[int], scientific_name: bool = False
    ) -> dict:
        if isinstance(tax_ids, int):
            tax_ids = [tax_ids]

        analysis_result = {
            "sample": {},
            "total_samples": 0,
            "rank_counts": defaultdict(
                lambda: defaultdict(lambda: {"count": 0, "tax_id": None})
            ),
        }

        for tax_id in tax_ids:
            # Get sample count and additional information
            sample_count = self._db.count(tax_id)
            tax_info = self._api[tax_id]
            scientific_name = tax_info.get("ScientificName", "Unknown")
            rank = tax_info.get("Rank", "Unknown")
            division = tax_info.get("Division", "Unknown")

            # Store sample information
            analysis_result["sample"][tax_id] = {
                "scientific_name": scientific_name,
                "rank": rank,
                "division": division,
                "count": sample_count,
            }
            analysis_result["total_samples"] += sample_count

            # Lineage information (ranks)
            lineage_ex = tax_info.get("LineageEx", [])

            # Count the number of elements for each parent rank
            for entry in lineage_ex:
                parent_rank = entry["Rank"]
                parent_scientific_name = entry["ScientificName"]
                parent_tax_id = int(entry["TaxId"])

                # Update the rank_counts structure
                analysis_result["rank_counts"][parent_rank][parent_scientific_name][
                    "count"
                ] += sample_count
                analysis_result["rank_counts"][parent_rank][parent_scientific_name][
                    "tax_id"
                ] = parent_tax_id

        # Convert defaultdict to a regular dict for the final result
        analysis_result["rank_counts"] = {
            k: dict(v) for k, v in analysis_result["rank_counts"].items()
        }

        return analysis_result
        rank_mapping = defaultdict(list)
        for tax_id in self._db.get_tax_ids():
            lineage_ex = self._api[tax_id].get("LineageEx", [])
            for entry in lineage_ex:
                if entry["Rank"] == rank:
                    rank_mapping[entry["TaxId"]].append(tax_id)
                    break

        return dict(rank_mapping)
