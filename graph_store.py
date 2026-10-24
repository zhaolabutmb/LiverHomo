"""Fast access to the PrePPI edge lists.

The edge lists are Excel workbooks with up to ~300,000 interactions; parsing
one takes up to two minutes. Each workbook is converted once to a compressed
NumPy file next to the other caches, which loads in well under a second.
"""

import os
import re
import zipfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Any, Iterator

import numpy as np

CACHE_VERSION = 1
_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_EDGE_COLUMNS = ("geneID1", "geneID2", "Score")


def _column_index(cell_ref: str) -> int:
    index = 0
    for char in cell_ref:
        if not char.isalpha():
            break
        index = index * 26 + ord(char.upper()) - 64
    return index - 1


def _xlsx_sheet_rows(path: str) -> Iterator[Iterator[list[Any]]]:
    """Yield the rows of each sheet. Streams the XML; about 3x faster than openpyxl."""
    with zipfile.ZipFile(path) as archive:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            for _, element in ET.iterparse(archive.open("xl/sharedStrings.xml")):
                if element.tag == _NS + "si":
                    shared.append("".join(text.text or "" for text in element.iter(_NS + "t")))
                    element.clear()
        sheet_names = [name for name in archive.namelist() if re.fullmatch(r"xl/worksheets/sheet\d+\.xml", name)]
        sheet_names.sort(key=lambda name: int(re.search(r"(\d+)\.xml$", name).group(1)))

        def rows(name: str) -> Iterator[list[Any]]:
            for _, element in ET.iterparse(archive.open(name)):
                if element.tag != _NS + "row":
                    continue
                values: list[Any] = []
                for cell in element.iter(_NS + "c"):
                    ref = cell.get("r")
                    if ref:
                        values.extend([None] * (_column_index(ref) - len(values)))
                    value = cell.find(_NS + "v")
                    cell_type = cell.get("t")
                    if value is None or value.text is None:
                        inline = cell.find(_NS + "is")
                        values.append("".join(text.text or "" for text in inline.iter(_NS + "t")) if inline is not None else None)
                    elif cell_type == "s":
                        values.append(shared[int(value.text)])
                    elif cell_type in ("str", "inlineStr", "e"):
                        values.append(value.text)
                    else:
                        values.append(float(value.text))
                element.clear()
                yield values

        for name in sheet_names:
            yield rows(name)


@dataclass
class GraphEdges:
    """Undirected interactions, strongest first. source/target index into genes."""

    genes: np.ndarray
    source: np.ndarray
    target: np.ndarray
    score: np.ndarray

    @property
    def edge_count(self) -> int:
        return int(self.score.size)

    def top(self, limit: int, seed_genes: list[str] | None = None) -> dict[str, Any]:
        """The `limit` highest-scoring edges, optionally only those touching seed_genes."""
        selected = np.arange(self.edge_count)
        if seed_genes:
            wanted = {gene.strip().upper() for gene in seed_genes if gene.strip()}
            seed_idx = np.flatnonzero(np.isin(np.char.upper(self.genes), list(wanted)))
            if seed_idx.size == 0:
                return {"nodes": [], "edges": [], "total_edges": 0}
            selected = np.flatnonzero(np.isin(self.source, seed_idx) | np.isin(self.target, seed_idx))
        total = int(selected.size)
        selected = selected[:limit]
        edges = [
            {"source": str(self.genes[s]), "target": str(self.genes[t]), "score": float(w)}
            for s, t, w in zip(self.source[selected], self.target[selected], self.score[selected])
        ]
        nodes = sorted({edge["source"] for edge in edges} | {edge["target"] for edge in edges})
        return {"nodes": nodes, "edges": edges, "total_edges": total}


def read_xlsx_edges(path: str) -> GraphEdges:
    """Combine all sheets; one edge per gene pair (highest score), no self-loops."""
    best: dict[tuple[str, str], float] = {}
    for rows in _xlsx_sheet_rows(path):
        headers = [str(value).strip() if value is not None else "" for value in next(rows, [])]
        if not all(name in headers for name in _EDGE_COLUMNS):
            continue
        columns = [headers.index(name) for name in _EDGE_COLUMNS]
        for row in rows:
            if len(row) <= max(columns):
                continue
            source, target, score = (row[index] for index in columns)
            if source is None or target is None or score is None:
                continue
            source, target = str(source).strip(), str(target).strip()
            if not source or not target or source == target:
                continue
            key = (source, target) if source < target else (target, source)
            score = float(score)
            if score > best.get(key, float("-inf")):
                best[key] = score

    genes = np.array(sorted({gene for pair in best for gene in pair}), dtype=str)
    position = {gene: index for index, gene in enumerate(genes)}
    pairs = sorted(best.items(), key=lambda item: (-item[1], item[0]))
    return GraphEdges(
        genes=genes,
        source=np.array([position[pair[0]] for pair, _ in pairs], dtype=np.int32),
        target=np.array([position[pair[1]] for pair, _ in pairs], dtype=np.int32),
        score=np.array([score for _, score in pairs], dtype=np.float32),
    )


def cache_path(xlsx_path: str, cache_dir: str) -> str:
    return os.path.join(cache_dir, os.path.splitext(os.path.basename(xlsx_path))[0] + ".npz")


def _source_stamp(xlsx_path: str) -> np.ndarray:
    stat = os.stat(xlsx_path)
    return np.array([CACHE_VERSION, stat.st_size, stat.st_mtime_ns], dtype=np.int64)


def load_edges(xlsx_path: str, cache_dir: str) -> GraphEdges:
    """Load from the cache, converting the workbook first if the cache is missing or stale."""
    target = cache_path(xlsx_path, cache_dir)
    stamp = _source_stamp(xlsx_path)
    if os.path.isfile(target):
        try:
            with np.load(target, allow_pickle=False) as data:
                if np.array_equal(data["stamp"], stamp):
                    return GraphEdges(data["genes"], data["source"], data["target"], data["score"])
        except (OSError, ValueError, KeyError):
            pass  # unreadable cache; rebuild it

    edges = read_xlsx_edges(xlsx_path)
    os.makedirs(cache_dir, exist_ok=True)
    temp_path = f"{target}.{os.getpid()}.tmp.npz"
    np.savez_compressed(temp_path, stamp=stamp, genes=edges.genes, source=edges.source, target=edges.target, score=edges.score)
    os.replace(temp_path, target)
    return edges
