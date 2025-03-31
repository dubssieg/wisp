from collections import defaultdict
import logging
import pickle
from typing import Any
import pandas as pd
from tqdm.auto import tqdm
from Bio import Entrez
from diskcache import Cache
from urllib.error import HTTPError
from pathlib import Path

LOG = logging.getLogger(__name__)

BACTERIA = 2  # root tax_id

RANKS = [
    "no rank",
    "cellular root",
    "domain",
    "superkingdom",
    "kingdom",
    "clade",
    "phylum",
    "class",
    "subclass",
    "order",
    "suborder",
    "family",
    "subfamily",
    "tribe",
    "genus",
    "subgenus",
    "species",
    "species group",
    "subspecies",
    "species subgroup",
    "strain",
    # ?
    "biotype",
    "pathogroup",
    "serogroup",
    "serotype",
]


class TaxDB:
    def __init__(
        self,
        cache_dir: str | Path,
        email: str,
        can_download: bool = True,
        preload: bool = False,
    ):
        LOG.debug(f"TaxDB({locals()})")
        self._can_download = can_download
        self._preload = preload
        self._cache_dir = Path(cache_dir).resolve()
        self._cache = Cache(self._cache_dir)
        self._tax_id_errors_ = set()
        Entrez.email = email
        if preload:
            self._pr_cache = {k: self._cache[k] for k in self._cache}

    def __getitem__(self, tax_id: int) -> dict | None:
        # preloaded
        if self._preload and tax_id in self._pr_cache:
            return self._pr_cache[tax_id][0]
        # in cache
        if record := self._cache.get(tax_id, None):
            return record[0]
        # not found, try do download
        if self._can_download:
            return self._get_api_data(tax_id)[0]
        # None
        return record

    def tax_ids(self) -> list:
        """ "List of all available tax_ids"""
        return list(self._cache)

    def populate_api_cache(self, metadata_csv_path: str | Path):
        """Run once."""
        errors = set()
        unique_tax_ids = list()
        extra_tax_ids = set()
        df = pd.read_csv(metadata_csv_path, sep="\t", usecols=["tax_id"])

        # direct tax_id
        for tax_id in tqdm(df["tax_id"].unique()):
            if self._get_api_data(tax_id):
                unique_tax_ids.append(tax_id)
            else:
                errors.add(tax_id)

        # extra tax_id
        for tax_id in tqdm(unique_tax_ids):
            if tax_id in self._cache:
                if tdata := self._cache[tax_id]:
                    extra_tax_ids.update(
                        [int(lin["TaxId"]) for lin in tdata[0].get("LineageEx", {})]
                    )
        for tax_id in tqdm(extra_tax_ids):
            self._get_api_data(tax_id)

        return unique_tax_ids, list(extra_tax_ids), list(errors)

    def clean_cache(self) -> list:
        """Clean the cache by removing entries with non-integer keys or empty values."""
        keys_to_delete = []

        for key in self._cache:
            value = self._cache.get(key)
            if not isinstance(key, int) or not value:
                keys_to_delete.append(key)

        for key in keys_to_delete:
            del self._cache[key]

        return keys_to_delete

    def export_db(self, pickle_path: str | Path | None = None):
        if pickle_path is None:
            pickle_path = Path(__file__).resolve().parent / "out/cache_dump.pkl"
            pickle_path.parent.mkdir(parents=True, exist_ok=True)
        pickle_path = Path(pickle_path)
        exp_cache = dict()
        for key in self._cache:
            exp_cache[key] = self._cache[key]

        with open(pickle_path, "wb") as f:
            pickle.dump(exp_cache, f)

    def import_db(self, pickle_path: str | Path | None = None):
        if pickle_path is None:
            pickle_path = Path(__file__).resolve().parent / "out/cache_dump.pkl"
        pickle_path = Path(pickle_path)
        with open(pickle_path, "rb") as f:
            exp_cache = pickle.load(f)
        for k, v in tqdm(exp_cache.items()):
            self._cache[k] = v

    def _get_api_data(self, tax_id: int | str) -> list | None:
        tax_id = self.clean_tax_id(tax_id)
        if not tax_id:
            return None

        # hit cache
        if tax_id in self._cache:
            return self._cache[tax_id]

        # API call + cache
        try:
            LOG.debug(f"Downloading taxonomy data for {tax_id}...")
            handle = Entrez.efetch(db="taxonomy", id=str(tax_id), retmode="xml")
            records = Entrez.read(handle)
            self._cache[tax_id] = records
            return records
        except HTTPError:
            LOG.error(f"Could not get API data for {tax_id}")
            return None

    def clean_tax_id(self, tax_id: Any) -> int:
        """Clean and normalize tax_id - find closest parent if needed"""

        # Clean the tax_id
        if not tax_id or pd.isna(tax_id):
            return BACTERIA

        if "," in str(tax_id):
            try:
                tids = [int(num) for num in tax_id.split(",")]
            except ValueError:
                return BACTERIA

        else:
            try:
                tids = [int(float(tax_id))]
            except ValueError:
                return BACTERIA

        if len(tids) == 1:
            return tids[0]

        # closest common parent tax_id
        lineage_map = defaultdict(lambda: (set(), None))

        for tid in tids:
            record = self[tid]
            lineage = [
                (entry["TaxId"], entry["Rank"]) for entry in record.get("LineageEx", [])
            ]
            # add self one tid is parent of an other tid
            lineage.append((record["TaxId"], record["Rank"]))
            if lineage:
                for ancestor_id, rank in lineage:
                    rank_index = RANKS.index(rank)
                    existing_rank_index = lineage_map[ancestor_id][1]
                    if existing_rank_index is None or rank_index > existing_rank_index:
                        lineage_map[ancestor_id] = (
                            lineage_map[ancestor_id][0],
                            rank_index,
                        )
                    lineage_map[ancestor_id][0].add(tid)

        common_ancestor = None
        highest_rank_index = -1
        for ancestor, (tids_set, rank_index) in lineage_map.items():
            if len(tids_set) == len(tids) and rank_index > highest_rank_index:
                common_ancestor = ancestor
                highest_rank_index = rank_index

        return int(common_ancestor) if common_ancestor else BACTERIA
