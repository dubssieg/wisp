from pathlib import Path
from typing import Literal

from diskcache import Cache
from functools import lru_cache

from reader import Reader

DB_TYPE = Literal["counter", "header", "index", "tmp"]


class Database:
    def __init__(self, kmer_size: int, dbs_path: str | Path, reader: Reader):
        self._dbs_path = Path(dbs_path).resolve()
        self._kmer_size = kmer_size
        self._reader = reader

    def push_file(self, file_path: str | Path):
        """Add content."""
        new_content = self._reader.process_file(
            file_path=file_path, kmer_size=self._kmer_size
        )
        print(new_content)

    @lru_cache(maxsize=None)
    def _get_sub_db(self, db_type: DB_TYPE, tax_id: int | None = None) -> Cache:
        """Get sub DB"""
        path = self._dbs_path / str(self._kmer_size) / db_type
        if tax_id:
            path /= str(tax_id)

        path.mkdir(parents=True, exist_ok=True)

        return Cache(path, timeout=None, size_limit=None)
