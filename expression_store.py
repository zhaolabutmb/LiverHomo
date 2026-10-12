import threading
from types import SimpleNamespace

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

# Rows scanned per step when extracting one gene from a CSR matrix. Bounds
# memory to a few hundred MB instead of loading the whole matrix.
COLUMN_SCAN_CHUNK = 32_000_000


class ExpressionStore:
    """Read-only access to the condensed gene-expression AnnData.

    AnnData's backed mode loads every layer into memory (~1.6 GB here), so
    this reads observation/variable tables eagerly and expression values one
    gene at a time directly from the HDF5 CSR matrices.
    """

    def __init__(self, path: str):
        self.file = h5py.File(path, "r")
        self.obs = read_elem(self.file["obs"])
        self.var = read_elem(self.file["var"])
        self.var_names = self.var.index
        self.n_obs, self.n_vars = len(self.obs), len(self.var)
        self.raw = None
        if "raw" in self.file:
            raw_var = read_elem(self.file["raw/var"])
            self.raw = SimpleNamespace(var=raw_var, var_names=raw_var.index)
        self.main_matrix = "layers/lognorm" if "layers/lognorm" in self.file else "X"
        self._indptr: dict[str, np.ndarray] = {}
        self._lock = threading.Lock()

    def _row_pointers(self, matrix: str) -> np.ndarray:
        with self._lock:
            if matrix not in self._indptr:
                self._indptr[matrix] = self.file[matrix]["indptr"][:]
            return self._indptr[matrix]

    def column(self, matrix: str, column_idx: int) -> np.ndarray:
        group = self.file[matrix]
        if group.attrs.get("encoding-type") != "csr_matrix":
            raise ValueError(f"{matrix} is not stored as a CSR matrix")
        indptr = self._row_pointers(matrix)
        indices, data = group["indices"], group["data"]
        values = np.zeros(indptr.size - 1, dtype=np.float32)
        for start in range(0, indices.shape[0], COLUMN_SCAN_CHUNK):
            hits = np.flatnonzero(indices[start:start + COLUMN_SCAN_CHUNK] == column_idx)
            if hits.size:
                positions = hits + start
                rows = np.searchsorted(indptr, positions, side="right") - 1
                values[rows] = data[positions[0]:positions[-1] + 1][positions - positions[0]]
        return values

    def close(self):
        self.file.close()


class MetadataStore:
    """Selected observation columns of an AnnData file, read into memory once."""

    def __init__(self, path: str, columns: tuple[str, ...]):
        with h5py.File(path, "r") as handle:
            obs = handle["obs"]
            present = [name for name in dict.fromkeys(columns) if name in obs]
            index = read_elem(obs[obs.attrs["_index"]])
            self.obs = pd.DataFrame({name: read_elem(obs[name]) for name in present}, index=pd.Index(index))
            self.total_obs_columns = len(obs.attrs.get("column-order", present))
            self.n_vars = handle["var/_index"].shape[0] if "var/_index" in handle else 0
        self.n_obs = len(self.obs)
