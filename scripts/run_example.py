"""Run the deterministic 28-point SNI-DPC example."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from sni_dpc import SNIDPC  # noqa: E402


def read_csv(path: Path) -> np.ndarray:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return np.asarray([[float(row["x"]), float(row["y"])] for row in rows])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "examples" / "data" / "worked_example_28.csv",
    )
    args = parser.parse_args()

    estimator = SNIDPC().fit(read_csv(args.input))
    summary = {
        "n_samples": int(estimator.labels_.size),
        "n_clusters": estimator.n_clusters_,
        "natural_neighbor_order": estimator.natural_neighbor_order_,
        "low_confidence_count": int(estimator.low_confidence_indices_.size),
        "prototype_sets": estimator.prototype_sets_,
    }
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

