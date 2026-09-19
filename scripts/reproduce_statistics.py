"""Recompute the historical v9.5 21-dataset summary and paired tests."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "results" / "synthetic21_six_method_results.csv"
METRICS = ("NMI", "ACC", "ARI")
METHODS = {
    "SNI-DPC": "MNG-DPC",
    "DPC-oracle-c": "DPC-oracle-c",
    "WANN-DPC": "WANNDPC",
    "DSNGCAP": "DSNGCAP",
    "Gravity+DSNGCAP": "Gravity_DSNGCAP",
    "KMeans++ (oracle-c)": "KNN",
}


def holm_adjust(p_values: list[float]) -> list[float]:
    order = sorted(range(len(p_values)), key=p_values.__getitem__)
    adjusted = [math.nan] * len(p_values)
    running_max = 0.0
    total = len(p_values)
    for rank, index in enumerate(order, start=1):
        running_max = max(running_max, min(1.0, (total - rank + 1) * p_values[index]))
        adjusted[index] = running_max
    return adjusted


def rank_sums(diff: np.ndarray) -> tuple[float, float, int]:
    nonzero = diff[np.abs(diff) > 1e-12]
    if not len(nonzero):
        return 0.0, 0.0, 0
    ranks = stats.rankdata(np.abs(nonzero), method="average")
    return float(ranks[nonzero > 0].sum()), float(ranks[nonzero < 0].sum()), len(nonzero)


def validate_input(frame: pd.DataFrame) -> None:
    if len(frame) != 21 or frame["dataset"].nunique() != 21:
        raise ValueError("The formal table must contain 21 unique datasets.")
    required = [f"{prefix}_{metric}" for prefix in METHODS.values() for metric in METRICS]
    missing = [column for column in required if column not in frame]
    if missing:
        raise ValueError(f"Missing result columns: {missing}")


def build_outputs(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    means = []
    tests = []
    raw_p = []
    for public_name, raw_prefix in METHODS.items():
        for metric in METRICS:
            values = frame[f"{raw_prefix}_{metric}"].astype(float)
            means.append({"method": public_name, "metric": metric, "mean": values.mean()})

    for metric in METRICS:
        baseline = frame[f"{METHODS['SNI-DPC']}_{metric}"].astype(float).to_numpy()
        for public_name, raw_prefix in METHODS.items():
            if public_name == "SNI-DPC":
                continue
            comparator = frame[f"{raw_prefix}_{metric}"].astype(float).to_numpy()
            diff = baseline - comparator
            wins = int((diff > 1e-12).sum())
            ties = int((np.abs(diff) <= 1e-12).sum())
            losses = int((diff < -1e-12).sum())
            w_plus, w_minus, n_nonzero = rank_sums(diff)
            if n_nonzero:
                result = stats.wilcoxon(
                    baseline,
                    comparator,
                    alternative="greater",
                    zero_method="wilcox",
                    method="auto",
                )
                p_value = float(result.pvalue)
            else:
                p_value = 1.0
            raw_p.append(p_value)
            tests.append(
                {
                    "metric": metric,
                    "comparator": public_name,
                    "wins": wins,
                    "ties": ties,
                    "losses": losses,
                    "mean_diff": float(diff.mean()),
                    "median_diff": float(np.median(diff)),
                    "W_plus": w_plus,
                    "W_minus": w_minus,
                    "n_nonzero": n_nonzero,
                    "p_one_sided": p_value,
                }
            )

    adjusted = holm_adjust(raw_p)
    for row, p_holm in zip(tests, adjusted):
        row["holm_p"] = p_holm
        row["significant"] = bool(p_holm < 0.05)
    return pd.DataFrame(means), pd.DataFrame(tests)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    args = parser.parse_args()

    frame = pd.read_csv(args.input)
    validate_input(frame)
    means, tests = build_outputs(frame)
    args.output.mkdir(parents=True, exist_ok=True)
    means.to_csv(args.output / "mean_metrics_21.csv", index=False)
    tests.to_csv(args.output / "wilcoxon_holm_21.csv", index=False)
    print(means.to_string(index=False))
    print(f"\nWrote {args.output / 'mean_metrics_21.csv'}")
    print(f"Wrote {args.output / 'wilcoxon_holm_21.csv'}")


if __name__ == "__main__":
    main()
