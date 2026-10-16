"""Reproduce seven-method score-level statistics; no clustering or label scoring.

Protocol seven-methods-2026-10-16: formal21 and real10 are separate inference
panels, with 24 joint-Holm hypotheses each. challenge8 is descriptive only.
The preserved six-method script supplies the previously reviewed mathematical
primitives; its inputs, entry point, and historical outputs remain unchanged.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd
import scipy
from scipy import stats

from reproduce_current_statistics import (
    canonical_sha256, exact_signed_rank_pvalue, hodges_lehmann_interval,
    holm_adjust,
)

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL_VERSION = "seven-methods-2026-10-16"
METHODS = ["SNI-DPC", "Gravity", "Torque Clustering", "LDCC", "Gauging-delta", "DPC-MFP", "WANN-DPC"]
METRICS = ["NMI", "ACC", "ARI", "Purity"]
EPS = 1e-12
ALPHA = 0.05
FAMILY_SIZE = 24


def load_panels(artificial: Path, real: Path) -> dict[str, pd.DataFrame]:
    manifest = json.loads((ROOT / "docs/benchmark-input-manifest.json").read_text(encoding="utf-8"))
    expected = {p: {r["dataset"] for r in manifest["datasets"] if r["panel"] == p}
                for p in ["formal21", "challenge8", "real10"]}
    frames = []
    for path, allowed in [(artificial, {"formal21", "challenge8"}), (real, {"real10"})]:
        frame = pd.read_csv(path, float_precision="round_trip")
        required = {"dataset", "method", "panel", "status", "true_c", "c_out", "reference_count_supplied", *METRICS}
        if required - set(frame.columns):
            raise ValueError("Missing columns: " + str(sorted(required - set(frame.columns))))
        if frame.duplicated(["dataset", "method"]).any():
            raise ValueError("Duplicate dataset/method rows")
        if set(frame.panel) != allowed or set(frame.method) != set(METHODS):
            raise ValueError("Seven-method input panels or methods differ from the protocol")
        for panel in allowed:
            part = frame[frame.panel == panel]
            if set(part.dataset) != expected[panel] or len(part) != len(METHODS) * len(expected[panel]):
                raise ValueError("Dataset coverage differs for " + panel)
            if any(set(g.method) != set(METHODS) for _, g in part.groupby("dataset")):
                raise ValueError("Incomplete method coverage")
        flags = frame.reference_count_supplied.map(lambda v: str(v).lower())
        if not flags.isin(["true", "false"]).all() or not flags.eq("true").equals(frame.method.eq("WANN-DPC")):
            raise ValueError("Reference-count access must identify WANN-DPC only")
        if not frame.status.isin(["ok", "timeout"]).all():
            raise ValueError("Unexpected execution state")
        for metric in METRICS:
            frame[metric] = pd.to_numeric(frame[metric], errors="raise")
            values = frame.loc[frame.status == "ok", metric].to_numpy(float)
            lower = -1.0 if metric == "ARI" else -EPS
            if not np.isfinite(values).all() or (values < lower).any() or (values > 1 + EPS).any():
                raise ValueError("Invalid successful metric: " + metric)
        if frame.loc[frame.status != "ok", [*METRICS, "c_out"]].notna().any().any():
            raise ValueError("Timeout quality/count fields must remain empty, never zero")
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True)
    missing = set(combined.loc[combined.status != "ok", ["dataset", "method"]].itertuples(index=False, name=None))
    if missing != {("page_blocks", "Gauging-delta")}:
        raise ValueError("The single registered timeout must be retained")
    return {panel: combined[combined.panel == panel].copy() for panel in expected}


def complete_pivot(frame: pd.DataFrame, metric: str) -> pd.DataFrame:
    part = frame[(frame.status == "ok") & frame[metric].notna()]
    return part.pivot(index="dataset", columns="method", values=metric).reindex(columns=METHODS).dropna(how="any")


def omnibus(panel, frame):
    records = []
    for metric in METRICS:
        table = complete_pivot(frame, metric)
        statistic, probability = stats.friedmanchisquare(*(table[m].to_numpy(float) for m in METHODS))
        records.append({"panel": panel, "metric": metric, "complete_datasets": len(table), "methods": len(METHODS),
            "friedman_chi_square": float(statistic), "friedman_p_asymptotic": float(probability),
            "kendall_w": float(statistic / (len(table) * (len(METHODS) - 1))),
            "sni_mean_rank_common_blocks": float(table.rank(axis=1, ascending=False, method="average")["SNI-DPC"].mean()),
            "calculation": "tie-corrected Friedman statistic with asymptotic chi-square probability"})
    return pd.DataFrame(records)


def pairwise(panel, frame):
    records, differences = [], []
    for metric in METRICS:
        part = frame[(frame.status == "ok") & frame[metric].notna()]
        table = part.pivot(index="dataset", columns="method", values=metric)
        for comparator in METHODS[1:]:
            matched = table.reindex(columns=["SNI-DPC", comparator]).dropna(how="any")
            diff = matched["SNI-DPC"].to_numpy(float) - matched[comparator].to_numpy(float)
            probability, positive, negative = exact_signed_rank_pvalue(diff)
            shift, low, high, critical = hodges_lehmann_interval(diff)
            denominator = positive + negative
            records.append({"panel": panel, "metric": metric, "reference": "SNI-DPC", "comparator": comparator,
                "reference_count_supplied_to_comparator": comparator == "WANN-DPC",
                "paired_datasets": len(diff), "nonzero_pairs": int((np.abs(diff) > EPS).sum()),
                "wins": int((diff > EPS).sum()), "ties": int((np.abs(diff) <= EPS).sum()), "losses": int((diff < -EPS).sum()),
                "mean_difference": float(diff.mean()), "median_difference": float(np.median(diff)),
                "hodges_lehmann_shift": shift, "hl_walsh_ci95_low": low, "hl_walsh_ci95_high": high,
                "walsh_critical_index_zero_based": critical,
                "rank_biserial": (positive - negative) / denominator if denominator else 0.0,
                "wilcoxon_w_plus": positive, "wilcoxon_w_minus": negative,
                "wilcoxon_p_two_sided_exact_signflip": probability,
                "zero_rule": "exclude abs(d)<=1e-12 from signed ranks; retain as ties",
                "tie_rule": "average ranks for tied absolute nonzero differences",
                "calculation": "exact conditional sign-permutation; no normal approximation",
                "ci_note": "nominal pointwise 95% Walsh-order interval; no-tie critical index; not multiplicity-adjusted; no bootstrap"})
            for dataset, value in zip(matched.index, diff):
                differences.append({"panel": panel, "dataset": dataset, "metric": metric, "reference": "SNI-DPC",
                    "comparator": comparator, "sni_score": float(matched.loc[dataset, "SNI-DPC"]),
                    "comparator_score": float(matched.loc[dataset, comparator]), "sni_minus_comparator": float(value),
                    "zero_under_protocol": bool(abs(value) <= EPS)})
    result = pd.DataFrame(records)
    if len(result) != FAMILY_SIZE:
        raise ValueError("The registered panel-wide correction family must contain 24 hypotheses")
    result["holm_family_size"] = FAMILY_SIZE
    result["holm_p_joint_24"] = holm_adjust(result.wilcoxon_p_two_sided_exact_signflip)
    result["holm_significant_joint_0_05"] = result.holm_p_joint_24 < ALPHA
    return result, pd.DataFrame(differences)


def descriptions(panels):
    records = []
    for panel, frame in panels.items():
        for metric in METRICS:
            common = complete_pivot(frame, metric)
            ranks = common.rank(axis=1, ascending=False, method="average")
            for method in METHODS:
                available = frame[(frame.method == method) & (frame.status == "ok")][metric].dropna()
                records.append({"panel": panel, "metric": metric, "method": method,
                    "reference_count_supplied": method == "WANN-DPC", "available_datasets": len(available),
                    "mean_score_available": float(available.mean()), "mean_score_common_blocks": float(common[method].mean()),
                    "complete_rank_blocks": len(common), "mean_rank_common_blocks": float(ranks[method].mean()),
                    "role": "descriptive_only" if panel == "challenge8" else "inferential_panel_description"})
    return pd.DataFrame(records)


def build_outputs(artificial: Path, real: Path):
    panels = load_panels(artificial, real)
    tests, diffs = zip(*(pairwise(p, panels[p]) for p in ["formal21", "real10"]))
    return {"friedman_seven.csv": pd.concat([omnibus(p, panels[p]) for p in ["formal21", "real10"]], ignore_index=True),
            "pairwise_seven.csv": pd.concat(tests, ignore_index=True),
            "descriptive_seven.csv": descriptions(panels),
            "paired_differences_seven.csv": pd.concat(diffs, ignore_index=True)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artificial", type=Path, default=ROOT / "results/artificial29_seven_method_results.csv")
    parser.add_argument("--real", type=Path, default=ROOT / "results/real10_seven_method_results.csv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs/seven_method_statistics")
    args = parser.parse_args()
    inputs = [args.artificial, args.real, ROOT / "docs/benchmark-input-manifest.json", ROOT / "docs/seven-method-evaluation-config.json"]
    names = ["friedman_seven.csv", "pairwise_seven.csv", "descriptive_seven.csv", "paired_differences_seven.csv", "protocol_seven.json"]
    if any((args.output_dir / name).resolve() in {p.resolve() for p in inputs} for name in names):
        parser.error("Output cannot overwrite an input")
    outputs = build_outputs(args.artificial, args.real)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, frame in outputs.items():
        frame.to_csv(args.output_dir / name, index=False, lineterminator="\n", float_format="%.17g")
    metadata = {"protocol_version": PROTOCOL_VERSION, "methods": METHODS, "metrics": METRICS,
        "inputs": {"artificial": {"file": args.artificial.name, "sha256_canonical_lf": canonical_sha256(args.artificial)},
                   "real": {"file": args.real.name, "sha256_canonical_lf": canonical_sha256(args.real)}},
        "script": {"file": "scripts/" + Path(__file__).name, "sha256_canonical_lf": canonical_sha256(Path(__file__))},
        "mathematical_dependency": {"file": "scripts/reproduce_current_statistics.py", "sha256_canonical_lf": canonical_sha256(ROOT / "scripts/reproduce_current_statistics.py")},
        "configuration": {"file": "docs/seven-method-evaluation-config.json", "sha256_canonical_lf": canonical_sha256(inputs[-1])},
        "dataset_identity": {"file": "docs/benchmark-input-manifest.json", "sha256_canonical_lf": canonical_sha256(inputs[-2])},
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
        "inferential_panels": ["formal21", "real10"], "descriptive_only": ["challenge8"],
        "friedman_complete_blocks": {"formal21": 21, "real10": 9}, "friedman_probability": "asymptotic chi-square",
        "wilcoxon": "two-sided exact conditional sign-permutation; zero tolerance 1e-12; average nonzero tie ranks",
        "holm_family": "24 hypotheses (4 metrics x 6 comparators), separately in each inferential panel",
        "interval": "nominal pointwise 95% Walsh-order with exact no-tie critical index; not multiplicity-adjusted or exact tie-adjusted; no bootstrap",
        "missingness": "retain the registered page_blocks/Gauging-delta timeout; no imputation",
        "scope": "score-level statistics only; no clustering, prediction-label scoring, or resource rerun",
        "outputs": {name: {"rows": len(frame), "sha256_canonical_lf": canonical_sha256(args.output_dir / name)} for name, frame in outputs.items()}}
    (args.output_dir / "protocol_seven.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"protocol": PROTOCOL_VERSION, "outputs": {k: len(v) for k, v in outputs.items()}, "clustering_calls": 0}))


if __name__ == "__main__":
    main()
