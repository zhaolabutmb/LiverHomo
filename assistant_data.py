"""Data lookups used by the website assistant's tools.

DEG files hold only significant genes (sig = up/down), one file per disease
and cell type. A gene missing from a file that exists means "tested, not
significant"; a missing file means that combination was not analysed.
"""
import json
import os
import re
import threading
from difflib import get_close_matches
from typing import Any

import pandas as pd

DISEASE_ALIASES = {
    "nash": "MASH", "nafld": "MASH", "masld": "MASH", "steatohepatitis": "MASH", "fatty liver": "MASH",
    "liver cancer": "HCC", "hepatocellular carcinoma": "HCC", "hepatoma": "HCC",
    "cholangiocarcinoma": "iCCA", "icc": "iCCA", "intrahepatic cholangiocarcinoma": "iCCA",
    "alcoholic liver disease": "ALD", "alcohol": "ALD", "alcoholic": "ALD",
    "primary biliary cholangitis": "PBC", "primary biliary cirrhosis": "PBC",
    "primary sclerosing cholangitis": "PSC",
    "portal vein tumor thrombus": "PVTT", "portal vein tumour thrombus": "PVTT",
    "metastatic lymph node": "MLN", "metastatic lymph nodes": "MLN", "lymph node metastasis": "MLN",
    "fetal liver": "Fetal", "foetal": "Fetal",
    "adjacent": "Adjacent Normal", "tumor adjacent": "Adjacent Normal", "adjacent normal": "Adjacent Normal",
    "healthy liver": "Healthy", "healthy control": "Healthy",
}

CELL_TYPE_GROUPS = {
    "t cell": ["CD4_CM", "CD4_EM", "CD4_Naive", "CD8_CM", "CD8_EM", "CD8_Naive", "CD8_TEMRA", "Treg", "MAIT", "gdT"],
    "cd4 t cell": ["CD4_CM", "CD4_EM", "CD4_Naive", "Treg"],
    "cd4": ["CD4_CM", "CD4_EM", "CD4_Naive", "Treg"],
    "cd8 t cell": ["CD8_CM", "CD8_EM", "CD8_Naive", "CD8_TEMRA"],
    "cd8": ["CD8_CM", "CD8_EM", "CD8_Naive", "CD8_TEMRA"],
    "nk cell": ["Circulating_NK_NKT", "Resident_NK"],
    "nk": ["Circulating_NK_NKT", "Resident_NK"],
    "dendritic cell": ["cDC1", "cDC2", "Mig.cDC", "pDC"],
    "dc": ["cDC1", "cDC2", "Mig.cDC", "pDC"],
    "myeloid": ["Macrophage", "Mono+mono_derived_cell", "Neutrophil", "Basophil", "cDC1", "cDC2", "Mig.cDC", "pDC"],
}

CELL_TYPE_ALIASES = {
    "hepatocyte": "Hepatocyte", "cholangiocyte": "Cholangiocyte", "biliary epithelial cell": "Cholangiocyte",
    "stellate": "Stellate_cell", "stellate cell": "Stellate_cell", "hsc": "Stellate_cell", "fibroblast": "Stellate_cell",
    "endothelial": "Endothelium", "endothelial cell": "Endothelium", "lsec": "Endothelium",
    "macrophage": "Macrophage", "kupffer": "Macrophage", "kupffer cell": "Macrophage",
    "monocyte": "Mono+mono_derived_cell", "mono": "Mono+mono_derived_cell",
    "b cell": "B_cell", "plasma cell": "Plasma", "neutrophil": "Neutrophil", "basophil": "Basophil",
    "regulatory t cell": "Treg", "treg": "Treg", "mait": "MAIT", "gamma delta t cell": "gdT", "gdt": "gdT",
    "plasmacytoid dendritic cell": "pDC", "migratory dc": "Mig.cDC", "migratory dendritic cell": "Mig.cDC",
}


def _key(text: str) -> str:
    text = re.sub(r"[_\-+./]", " ", str(text).lower())
    text = re.sub(r"\s+", " ", text).strip()
    return re.sub(r"s$", "", text) if len(text) > 3 else text


class NameResolver:
    def __init__(self, diseases: dict[str, str], cell_types: list[str]):
        # diseases: canonical short name -> full display name
        self.diseases = diseases
        self.cell_types = cell_types
        self._disease_keys = {_key(name): name for name in diseases}
        self._disease_keys.update({_key(full): name for name, full in diseases.items()})
        self._disease_keys.update({_key(alias): name for alias, name in DISEASE_ALIASES.items() if name in diseases})
        self._cell_keys = {_key(name): [name] for name in cell_types}
        self._cell_keys.update({_key(alias): [name] for alias, name in CELL_TYPE_ALIASES.items() if name in cell_types})
        self._cell_keys.update({_key(alias): [n for n in names if n in cell_types] for alias, names in CELL_TYPE_GROUPS.items()})

    def disease(self, text: str | None) -> str | None:
        if not text or not str(text).strip():
            return None
        key = _key(text)
        match = self._disease_keys.get(key) or self._disease_keys.get(_key(re.sub(r"\(.*?\)", "", text)))
        if match is None:
            inner = re.search(r"\(([^)]+)\)", text)
            match = self._disease_keys.get(_key(inner.group(1))) if inner else None
        if match is None:
            raise LookupError(f"Unknown disease '{text}'. Available: {', '.join(self.diseases)}")
        return match

    def cell_types_for(self, text: str | None) -> list[str]:
        if not text or not str(text).strip():
            return []
        resolved: list[str] = []
        for part in re.split(r",|;|\band\b", str(text)):
            if not part.strip():
                continue
            key = _key(part)
            names = self._cell_keys.get(key) or self._cell_keys.get(re.sub(r" cell$", "", key))
            if not names:
                close = get_close_matches(key, list(self._cell_keys), n=3, cutoff=0.75)
                if len(close) == 1:
                    names = self._cell_keys[close[0]]
            if not names:
                raise LookupError(f"Unknown cell type '{part.strip()}'. Available: {', '.join(self.cell_types)}")
            resolved.extend(name for name in names if name not in resolved)
        return resolved


class DegIndex:
    COLUMNS = ["gene", "de_coef", "de_pval", "sig", "protein"]

    def __init__(self, deg_root: str):
        self.deg_root = deg_root
        self._table: pd.DataFrame | None = None
        self._tested: set[tuple[str, str]] = set()
        self._lock = threading.Lock()

    @property
    def table(self) -> pd.DataFrame:
        with self._lock:
            if self._table is None:
                self._table = self._load()
            return self._table

    @property
    def tested(self) -> set[tuple[str, str]]:
        self.table
        return self._tested

    def _load(self) -> pd.DataFrame:
        frames = []
        for disease in sorted(os.listdir(self.deg_root)):
            folder = os.path.join(self.deg_root, disease)
            if not os.path.isdir(folder):
                continue
            prefix = f"DEG_{disease}_"
            for filename in sorted(os.listdir(folder)):
                if not (filename.startswith(prefix) and filename.endswith(".csv")):
                    continue
                cell_type = filename[len(prefix):-4]
                self._tested.add((disease, cell_type))
                frame = pd.read_csv(os.path.join(folder, filename), sep="\t", usecols=self.COLUMNS)
                frame = frame[frame["sig"].isin(["up", "down"])]
                if frame.empty:
                    continue
                frame["disease"] = disease
                frame["cell_type"] = cell_type
                frames.append(frame)
        table = pd.concat(frames, ignore_index=True)
        for column in ("disease", "cell_type", "sig", "protein"):
            table[column] = table[column].astype("category")
        table["gene_upper"] = table["gene"].str.upper()
        return table

    @staticmethod
    def _rows(frame: pd.DataFrame, limit: int) -> list[dict[str, Any]]:
        return [
            {
                "gene": row.gene, "disease": row.disease, "cell_type": row.cell_type, "direction": row.sig,
                "de_coef": round(float(row.de_coef), 3),
                "de_pval": float(f"{row.de_pval:.3g}") if row.de_pval > 0 else "< 1e-300",
            }
            for row in frame.head(limit).itertuples(index=False)
        ]

    def gene(self, gene: str, diseases: list[str], cell_types: list[str], limit: int = 40) -> dict[str, Any]:
        table = self.table
        hits = table[table["gene_upper"] == gene.strip().upper()]
        if diseases:
            hits = hits[hits["disease"].isin(diseases)]
        if cell_types:
            hits = hits[hits["cell_type"].isin(cell_types)]
        hits = hits.sort_values("de_pval")
        tested = sorted(
            f"{d} / {c}" for d, c in self.tested
            if (not diseases or d in diseases) and (not cell_types or c in cell_types)
        )
        significant = {f"{d} / {c}" for d, c in zip(hits["disease"], hits["cell_type"])}
        result: dict[str, Any] = {
            "gene": gene.strip().upper() if hits.empty else hits["gene"].iloc[0],
            "protein": None if hits.empty else hits["protein"].iloc[0],
            "significant_hits": int(len(hits)),
            "up": int((hits["sig"] == "up").sum()),
            "down": int((hits["sig"] == "down").sum()),
            "rows": self._rows(hits, limit),
            "tested_combinations": len(tested),
        }
        not_significant = [combo for combo in tested if combo not in significant]
        if len(not_significant) <= 30:
            result["tested_not_significant"] = not_significant
        if not tested:
            result["note"] = "No DEG analysis exists for that disease/cell type selection."
        return result

    def top(self, disease: str, cell_types: list[str], direction: str | None, limit: int, sort: str) -> dict[str, Any]:
        table = self.table
        hits = table[table["disease"] == disease]
        if cell_types:
            hits = hits[hits["cell_type"].isin(cell_types)]
        if direction in ("up", "down"):
            hits = hits[hits["sig"] == direction]
        if sort == "effect":
            hits = hits.reindex(hits["de_coef"].abs().sort_values(ascending=False).index)
        else:
            hits = hits.sort_values(["de_pval", "de_coef"], ascending=[True, direction == "down"])
        missing = [c for c in cell_types if (disease, c) not in self.tested]
        result = {
            "disease": disease, "cell_types": cell_types or "all", "direction": direction or "both",
            "sorted_by": "absolute de_coef" if sort == "effect" else "de_pval",
            "total_significant": int(len(hits)), "rows": self._rows(hits, limit),
        }
        if missing:
            result["not_analysed"] = missing
        return result


def load_composition(path: str) -> dict[str, dict[str, int]]:
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def short_disease_name(full_name: str) -> str:
    inner = re.search(r"\(([^)]+)\)\s*$", full_name)
    return inner.group(1) if inner else full_name
