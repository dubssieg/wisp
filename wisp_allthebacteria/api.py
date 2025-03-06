import traceback
import pandas as pd
import tqdm
from Bio import Entrez
from diskcache import Cache
from urllib.error import HTTPError
from pathlib import Path

API_CACHE_DIR = "/data/microtaxo/apicache"
EMAIL = "cyrille.leroux@irisa.fr"


class API:
    def __init__(self, can_download: bool = True):
        self._can_download = can_download
        self._api_cache_dir = API_CACHE_DIR
        self._api_cache = Cache(self._api_cache_dir)
        self._tax_id_errors_ = set()
        Entrez.email = EMAIL

    def __getitem__(self, tax_id: int) -> dict:
        if record := self._api_cache.get(tax_id, None):
            return record[0]
        if self._can_download:
            return self._get_api_data(tax_id)[0]
        return record

    def tax_ids(self) -> list:
        """ "List of all available tax_ids"""
        return list(self._api_cache)

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
            if tax_id in self._api_cache:
                if tdata := self._api_cache[tax_id]:
                    extra_tax_ids.update(
                        [int(lin["TaxId"]) for lin in tdata[0].get("LineageEx", {})]
                    )
        for tax_id in tqdm(extra_tax_ids):
            self._get_api_data(tax_id)

        return unique_tax_ids, list(extra_tax_ids), list(errors)

    @staticmethod
    def clean_tax_id(tax_id) -> list:
        if pd.isna(tax_id):
            return []

        if "," in str(tax_id):
            try:
                return [int(num) for num in tax_id.split(",")]
            except ValueError:
                return []

        try:
            return [int(float(tax_id))]
        except ValueError:
            return []

    def clean_cache(self) -> list:
        """Clean the cache by removing entries with non-integer keys or empty values."""
        keys_to_delete = []

        for key in self._api_cache:
            value = self._api_cache.get(key)
            if not isinstance(key, int) or not value:
                keys_to_delete.append(key)

        for key in keys_to_delete:
            del self._api_cache[key]

        return keys_to_delete

    def _get_api_data(self, tax_id):
        tax_id = self.clean_tax_id(tax_id)
        if not tax_id:
            return None
        tax_id = tax_id[0]

        # hit cache
        if tax_id in self._api_cache:
            return self._api_cache[tax_id]

        # API call + cache
        try:
            handle = Entrez.efetch(db="taxonomy", id=str(tax_id), retmode="xml")
            records = Entrez.read(handle)
            self._api_cache[tax_id] = records
            return records
        except HTTPError:
            traceback.print_exc()
            return None
