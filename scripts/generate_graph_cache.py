"""Convert the PPI edge-list workbooks in GRAPH_ROOT to fast .npz files in GRAPH_CACHE_DIR.

The app also converts a workbook on first use, but that takes up to two
minutes per network, so run this after adding or replacing edge lists.
"""

import os
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from env_config import load_env  # noqa: E402
from graph_store import load_edges  # noqa: E402

load_env()


def main() -> None:
    graph_root = os.environ["GRAPH_ROOT"]
    cache_dir = os.environ.get("GRAPH_CACHE_DIR") or os.path.join(
        os.path.dirname(os.environ["GENE_VECTOR_CACHE_DIR"]), "ppi_edges"
    )
    paths = sorted(Path(graph_root).glob("*.xlsx"))
    print(f"{len(paths)} edge lists -> {cache_dir}", flush=True)
    started = time.time()
    for number, path in enumerate(paths, start=1):
        step = time.time()
        edges = load_edges(str(path), cache_dir)
        print(f"[{number}/{len(paths)}] {path.name}: {edges.edge_count} edges, {time.time() - step:.1f} s", flush=True)
    print(f"Done in {(time.time() - started) / 60:.1f} min")


if __name__ == "__main__":
    main()
