import logging
from pathlib import Path
from typing import Generator
import xgboost as xgb


LOG = logging.getLogger(__name__)
DMAT_SUFFIX = ".dmat"


class DMatStore:
    def __init__(self, path: str | Path, signature: str | None = None):
        self._path = Path(path).resolve()
        self._signature = str(signature)

    def has_dmat(self, name: str) -> bool:
        return self._dmat_path(name).exists()

    def dmat_list(self) -> list[str]:
        return [p.name for p in self._storage_path().rglob(f"*{DMAT_SUFFIX}")]

    def clear(self) -> None:
        for dmat_name in self.dmat_list():
            LOG.debug(f"deleting {dmat_name}...")
            dmat_path = self._dmat_path(dmat_name)
            dmat_path.unlink()

    def dmat_generator(self, names: list[str]) -> Generator[xgb.DMatrix, None, None]:
        for name in names:
            yield self._deserialize_dmatrix(name)

    def _storage_path(self) -> Path:
        path = self._path
        if self._signature:
            path /= self._signature
        return path

    def _dmat_path(self, name: str | int) -> Path:
        return self._storage_path() / f"{name}{DMAT_SUFFIX}"

    def _serialize_dmatrix(self, dmat: xgb.DMatrix, name: str) -> None:
        dmat_path = self._dmat_path(name)
        dmat_path.parent.mkdir(parents=True, exist_ok=True)
        LOG.debug(f"Storing DMatrix to {dmat_path}")
        dmat.save_binary(dmat_path)

    def _deserialize_dmatrix(self, name: str) -> xgb.DMatrix:
        dmat_path = self._dmat_path(name)
        if not dmat_path.exists():
            raise FileNotFoundError(f"{dmat_path} DMatrix not found.")
        return xgb.DMatrix(dmat_path)

    def __getitem__(self, name: str) -> xgb.DMatrix:
        return self._deserialize_dmatrix(name)

    def __setitem__(self, name: str, dmat: xgb.DMatrix):
        self._serialize_dmatrix(name=name, dmat=dmat)
