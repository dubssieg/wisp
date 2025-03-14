import pandas as pd
from api import API
from pathlib import Path


class Metadata:
    def __init__(self, csv_path: str | Path, api: API, start_loaded: bool = False):
        self._csv_path = Path(csv_path)
        self._md = None
        self._api = api
        if start_loaded:
            self.md

    @property
    def md(self):
        if self._md is None:
            self._md = self._load()
        return self._md

    def __getitem__(self, seq_id: str) -> dict:
        """Example: SAMD00013333 or SAMD00013333.contig0000"""
        if ".contig" in seq_id:
            seq_id = seq_id.split(".contig")[0]
        data = self.md[seq_id]
        if not data:
            return []

        if not isinstance(data[0], dict):
            data = [self._api[d] for d in data]
            self._md[seq_id] = data
        return data

    def _load(self) -> dict:
        """Read DataFrame once then free memory (big file)."""
        columns = ["sample_accession", "tax_id"]
        dtype_specification = {col: str for col in columns}
        df = pd.read_csv(
            self._csv_path, sep="\t", usecols=columns, dtype=dtype_specification
        )
        return {
            seq_id: API.clean_tax_id(tax_id)
            for seq_id, tax_id in zip(df[columns[0]], df[columns[1]])
        }
