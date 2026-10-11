import json
import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.chdir(PROJECT_ROOT)
sys.path.insert(0, str(PROJECT_ROOT))

from app import (  # noqa: E402
    get_cell_interaction_catalog,
    get_cell_interaction_cell_types,
    get_cell_type_options,
    get_deg_summary,
    get_disease_options,
    get_ego_catalog,
    get_graph_summary,
    get_volcano_catalog,
    resolve_cell_type_column,
)


def main() -> None:
    catalog = {
        "metadata": {
            "diseases": get_disease_options(),
            "cell_types": get_cell_type_options(resolve_cell_type_column()),
        },
        "deg": get_deg_summary(),
        "volcano_diseases": [item["disease"] for item in get_volcano_catalog()],
        "graph": get_graph_summary(),
        "ego": get_ego_catalog(),
        "cell_interactions": {
            "diseases": [item["disease"] for item in get_cell_interaction_catalog("heatmap")],
            "cell_types": get_cell_interaction_cell_types(),
        },
    }
    output_path = PROJECT_ROOT / "static" / "filter_catalog.json"
    output_path.write_text(json.dumps(catalog, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Generated selector catalog: {output_path}")


if __name__ == "__main__":
    main()
