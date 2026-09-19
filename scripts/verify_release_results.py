"""Validate the published 29-artificial + 10-real six-method result panels.

The tables are frozen experiment records. This script verifies coverage and
missing-value semantics; it does not rerun any algorithm or recompute metrics.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METHODS = {
    "SNI-DPC",
    "Gravity",
    "Torque Clustering",
    "LDCC",
    "Gauging-delta",
    "DPC-MFP",
}
METRICS = ("NMI", "ACC", "ARI", "Purity")
PANELS = (
    (ROOT / "results" / "artificial29_six_method_results.csv", 29),
    (ROOT / "results" / "real10_six_method_results.csv", 10),
)


def check_panel(path: Path, expected_datasets: int) -> tuple[int, int, int]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    grouped: dict[str, dict[str, dict[str, str]]] = {}
    timeout_count = 0
    for row in rows:
        dataset, method, status = row["dataset"], row["method"], row["status"]
        if not dataset or method not in METHODS:
            raise ValueError(f"Invalid dataset/method: {dataset!r}/{method!r}")
        group = grouped.setdefault(dataset, {})
        if method in group:
            raise ValueError(f"Duplicate result: {dataset}/{method}")
        group[method] = row
        if status == "ok":
            for metric in METRICS:
                value = float(row[metric])
                lower = -1.0 if metric == "ARI" else -1e-12
                if not math.isfinite(value) or not lower <= value <= 1.0 + 1e-12:
                    raise ValueError(f"Invalid {metric}: {dataset}/{method}")
            if int(float(row["c_out"])) < 1:
                raise ValueError(f"Invalid cluster count: {dataset}/{method}")
        elif status == "timeout":
            timeout_count += 1
            if any(row[metric].strip() for metric in METRICS) or row["c_out"].strip():
                raise ValueError(f"Timeout must not have quality/count: {dataset}/{method}")
        else:
            raise ValueError(f"Unrecognized status: {dataset}/{method}/{status}")
    if len(grouped) != expected_datasets:
        raise ValueError(f"Expected {expected_datasets} datasets in {path}, found {len(grouped)}")
    for dataset, group in grouped.items():
        if set(group) != METHODS:
            raise ValueError(f"Incomplete methods on {dataset}: {METHODS - set(group)}")
        if len({row["true_c"] for row in group.values()}) != 1:
            raise ValueError(f"Inconsistent reference count on {dataset}")
    if len(rows) != expected_datasets * len(METHODS):
        raise ValueError(f"Unexpected number of records in {path}")
    return len(grouped), len(rows), timeout_count


def main() -> None:
    for path, expected in PANELS:
        datasets, rows, timeouts = check_panel(path, expected)
        print(f"{path.name}: {datasets} datasets, {rows} records, {timeouts} timeout(s)")


if __name__ == "__main__":
    main()
