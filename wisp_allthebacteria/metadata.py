import logging
from typing import Any
import pandas as pd
from taxdb import TaxDB
from pathlib import Path

LOG = logging.getLogger(__name__)


class Metadata:
    def __init__(
        self,
        csv_path: str | Path,
        taxdb: TaxDB,
        start_loaded: bool = False,
        keep_data: bool = False,
    ):
        LOG.debug(f"Metadata({locals()})")
        self._csv_path = Path(csv_path)
        self._md = None
        self._atb_md = None  # AllTheBacteria MD
        self._taxdb = taxdb
        self._keep_data = keep_data
        self._dm = None
        if start_loaded:
            self.md

    @property
    def md(self) -> dict:
        if self._md is None:
            self._md = self._load()
        return self._md

    @property
    def atb_md(self) -> dict:
        if self._atb_md is None:
            self._atb_md = self._atb_md_load()
        return self._atb_md

    def clean_tax_id(self, tax_id: Any) -> int:
        """Expose TaxDB clean_tax_id"""
        return self._taxdb.clean_tax_id(tax_id)

    def get_atb_md(self, tax_id: int) -> list:
        matching_rows = self.atb_md.loc[self.atb_md["tax_id"] == tax_id]
        if matching_rows.empty:
            return []
        return matching_rows.to_dict(orient="records")

    def __getitem__(self, seq_id: str) -> dict:
        """Example: SAMD00013333 or SAMD00013333.contig0000
        Unusual exemple: SAMEA3924086.NODE_1_length_606878_cov_39.021295_pilon"""
        if "." in seq_id:
            seq_id = seq_id.split(".")[0]
        data = self.md[seq_id]
        if not data:
            return {}

        if self._keep_data:
            if not isinstance(data, dict):
                data = self._taxdb[data]
                self._md[seq_id] = data
        else:
            data = self._taxdb[data]
        return data

    def _clean_seq_id(self, seq_id: str) -> str:
        if "." in seq_id:
            seq_id = seq_id.split(".")[0]
        return seq_id

    def _load(self) -> dict:
        """Read DataFrame once then free memory (big file)."""
        columns = ["sample_accession", "tax_id"]
        LOG.debug(f"Reading {self._csv_path}...")
        dtype_specification = {col: str for col in columns}
        df = pd.read_csv(
            self._csv_path, sep="\t", usecols=columns, dtype=dtype_specification
        )
        return {
            seq_id: self.clean_tax_id(tax_id)
            for seq_id, tax_id in zip(df[columns[0]], df[columns[1]])
        }

    def _atb_md_load(self) -> dict:
        """Read DataFrame and keep it."""
        LOG.debug(f"Reading {self._csv_path} (FULL)...")
        df = pd.read_csv(self._csv_path, sep="\t")
        LOG.debug("Cleaning tax_ids...")
        df["tax_id_cleaned"] = df["tax_id"].apply(self.clean_tax_id)
        return df
