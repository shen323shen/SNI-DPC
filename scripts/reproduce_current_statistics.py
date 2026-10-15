"""Reproduce current six-method statistics from published result tables only.

Protocol: current-six-methods-2026-10-15. The original 21 artificial datasets
and ten real datasets form separate inferential panels. The eight artificial
extensions are descriptive only. No clustering, metric recomputation from
labels, bootstrap, or resource measurement is performed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL_VERSION = "current-six-methods-2026-10-15"
METHODS = ["SNI-DPC", "Gravity", "Torque Clustering", "LDCC", "Gauging-delta", "DPC-MFP"]
METRICS = ["NMI", "ACC", "ARI", "Purity"]
EPS = 1e-12
ALPHA = 0.05


def canonical_sha256(path: Path) -> str:
    """Hash UTF-8 text with CRLF converted to LF, matching release checksums."""
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_panels(artificial: Path, real: Path) -> dict[str, pd.DataFrame]:
    manifest = json.loads((ROOT / "docs" / "benchmark-input-manifest.json").read_text(encoding="utf-8"))
    expected = {
        panel: {row["dataset"] for row in manifest["datasets"] if row["panel"] == panel}
        for panel in ("formal21", "challenge8", "real10")
    }
    frames = []
    for path, panels in ((artificial, ("formal21", "challenge8")), (real, ("real10",))):
        frame = pd.read_csv(path)
        required = {"dataset", "method", "status", "panel", "true_c", "c_out", *METRICS}
        if required.difference(frame.columns):
            raise ValueError(f"{path.name}: missing columns {sorted(required.difference(frame.columns))}")
        if frame.duplicated(["dataset", "method"]).any():
            raise ValueError(f"{path.name}: duplicate dataset/method records")
        if panels == ("real10",):
            # The frozen public table retains the collection's original panel name.
            if set(frame["panel"]) != {"real_top30"}:
                raise ValueError(f"{path.name}: expected frozen real_top30 source panel")
            frame["panel"] = "real10"
        if set(frame["panel"]) != set(panels) or set(frame["method"]) != set(METHODS):
            raise ValueError(f"{path.name}: unexpected panels or methods (historical input is incompatible)")
        for panel in panels:
            subset = frame[frame["panel"] == panel]
            if set(subset["dataset"]) != expected[panel] or len(subset) != 6 * len(expected[panel]):
                raise ValueError(f"{path.name}: dataset coverage mismatch for {panel}")
            if any(set(group["method"]) != set(METHODS) for _, group in subset.groupby("dataset")):
                raise ValueError(f"{path.name}: incomplete method coverage for {panel}")
        if not frame["status"].isin(["ok", "timeout"]).all():
            raise ValueError(f"{path.name}: unknown non-success state")
        for metric in METRICS:
            frame[metric] = pd.to_numeric(frame[metric], errors="raise")
            valid = frame.loc[frame["status"] == "ok", metric].to_numpy(float)
            lower = -1.0 if metric == "ARI" else -EPS
            if not np.isfinite(valid).all() or (valid < lower).any() or (valid > 1 + EPS).any():
                raise ValueError(f"{path.name}: invalid successful {metric}")
        if frame.loc[frame["status"] != "ok", [*METRICS, "c_out"]].notna().any().any():
            raise ValueError(f"{path.name}: non-success quality/count must remain missing, not zero")
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True)
    timeouts = set(combined.loc[combined["status"] == "timeout", ["dataset", "method"]].itertuples(index=False, name=None))
    if timeouts != {("page_blocks", "Gauging-delta")}:
        raise ValueError("Current frozen tables require exactly the recorded page_blocks/Gauging-delta timeout")
    return {panel: combined[combined["panel"] == panel].copy() for panel in expected}


def complete_pivot(frame: pd.DataFrame, metric: str) -> pd.DataFrame:
    valid = frame[(frame["status"] == "ok") & frame[metric].notna()]
    return valid.pivot(index="dataset", columns="method", values=metric).reindex(columns=METHODS).dropna(how="any")


def exact_signed_rank_pvalue(differences: np.ndarray) -> tuple[float, float, float]:
    values = np.asarray(differences, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Signed-rank differences must be finite")
    nonzero = values[np.abs(values) > EPS]
    if len(nonzero) == 0:
        return 1.0, 0.0, 0.0
    ranks = stats.rankdata(np.abs(nonzero), method="average")
    scaled = np.rint(2.0 * ranks).astype(int)
    observed = int(np.rint(2.0 * ranks[nonzero > 0].sum()))
    counts = [0] * (int(scaled.sum()) + 1)
    counts[0] = 1
    reachable = 0
    for rank in scaled:
        rank = int(rank)
        for subtotal in range(reachable, -1, -1):
            counts[subtotal + rank] += counts[subtotal]
        reachable += rank
    assignments = 2 ** len(scaled)
    lower = sum(counts[: observed + 1]) / assignments
    upper = sum(counts[observed:]) / assignments
    return float(min(1.0, 2 * min(lower, upper))), float(ranks[nonzero > 0].sum()), float(ranks[nonzero < 0].sum())


def holm_adjust(values) -> np.ndarray:
    raw = np.asarray(values, dtype=float)
    if not np.isfinite(raw).all() or (raw < 0).any() or (raw > 1).any():
        raise ValueError("Holm inputs must be finite probabilities")
    adjusted = np.empty_like(raw)
    running = 0.0
    for position, index in enumerate(np.argsort(raw)):
        running = max(running, min(1.0, (len(raw) - position) * raw[index]))
        adjusted[index] = running
    return adjusted


def signed_rank_critical_value(n: int, tail_probability: float) -> int:
    counts = [0] * (n * (n + 1) // 2 + 1)
    counts[0] = 1
    reachable = 0
    for rank in range(1, n + 1):
        for subtotal in range(reachable, -1, -1):
            counts[subtotal + rank] += counts[subtotal]
        reachable += rank
    cumulative, critical = 0, -1
    for value, count in enumerate(counts):
        cumulative += count
        if cumulative / 2**n > tail_probability:
            break
        critical = value
    return critical


def hodges_lehmann_interval(differences: np.ndarray) -> tuple[float, float, float, int]:
    values = np.asarray(differences, dtype=float)
    if not len(values) or not np.isfinite(values).all():
        raise ValueError("Walsh interval requires at least one finite pair")
    # Preserve the reviewed protocol: all pairs enter Walsh averages, including zeros.
    walsh = np.sort([(values[i] + values[j]) / 2 for i in range(len(values)) for j in range(i, len(values))])
    critical = signed_rank_critical_value(len(values), ALPHA / 2)
    estimate = float(np.median(walsh))
    if critical < 0:
        return estimate, float("nan"), float("nan"), critical
    return estimate, float(walsh[critical]), float(walsh[len(walsh) - critical - 1]), critical


def friedman_rows(panel: str, frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in METRICS:
        pivot = complete_pivot(frame, metric)
        statistic, p_value = stats.friedmanchisquare(*(pivot[method].to_numpy(float) for method in METHODS))
        n = len(pivot)
        rows.append(dict(panel=panel, metric=metric, complete_datasets=n, methods=6,
                         friedman_chi_square=float(statistic), friedman_p=float(p_value),
                         kendall_w=float(statistic / (n * 5)),
                         inferential_role="primary_artificial" if panel == "formal21" else "separate_real_data"))
    return pd.DataFrame(rows)


def pairwise_rows(panel: str, frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in METRICS:
        valid = frame[(frame["status"] == "ok") & frame[metric].notna()]
        pivot = valid.pivot(index="dataset", columns="method", values=metric)
        for comparator in METHODS[1:]:
            pair = pivot.reindex(columns=["SNI-DPC", comparator]).dropna(how="any")
            diff = pair["SNI-DPC"].to_numpy(float) - pair[comparator].to_numpy(float)
            p_value, w_plus, w_minus = exact_signed_rank_pvalue(diff)
            hl, low, high, critical = hodges_lehmann_interval(diff)
            denominator = w_plus + w_minus
            rows.append(dict(
                panel=panel, metric=metric, reference="SNI-DPC", comparator=comparator,
                paired_datasets=len(pair), nonzero_pairs=int((np.abs(diff) > EPS).sum()),
                wins=int((diff > EPS).sum()), ties=int((np.abs(diff) <= EPS).sum()), losses=int((diff < -EPS).sum()),
                mean_difference=float(diff.mean()), median_difference=float(np.median(diff)),
                hodges_lehmann_shift=hl, hl_walsh_ci95_low=low, hl_walsh_ci95_high=high,
                walsh_critical_index_zero_based=critical,
                rank_biserial=(w_plus - w_minus) / denominator if denominator else 0.0,
                wilcoxon_w_plus=w_plus, wilcoxon_w_minus=w_minus,
                wilcoxon_p_two_sided_exact_signflip=p_value,
                zero_rule="exclude abs(d)<=1e-12 from signed ranks; retain as ties",
                tie_rule="average ranks for tied absolute nonzero differences",
                calculation="exact conditional sign-permutation; not asymptotic",
                ci_note="Walsh-order signed-rank 95% interval; no bootstrap; conservative with ties",
                inferential_role="primary_artificial" if panel == "formal21" else "separate_real_data"))
    result = pd.DataFrame(rows)
    result["holm_p_joint_20"] = holm_adjust(result["wilcoxon_p_two_sided_exact_signflip"])
    result["holm_significant_joint_0_05"] = result["holm_p_joint_20"] < ALPHA
    result["holm_p_within_metric_5_audit_only"] = result.groupby("metric", sort=False)[
        "wilcoxon_p_two_sided_exact_signflip"].transform(holm_adjust)
    return result


def descriptive_rows(panels: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for panel, frame in panels.items():
        for metric in METRICS:
            complete = complete_pivot(frame, metric)
            ranks = complete.rank(axis=1, ascending=False, method="average")
            for method in METHODS:
                available = frame[(frame["method"] == method) & (frame["status"] == "ok")][metric].dropna()
                rows.append(dict(panel=panel, metric=metric, method=method, available_datasets=len(available),
                                 mean_score=float(available.mean()), complete_rank_blocks=len(complete),
                                 mean_rank_common_blocks=float(ranks[method].mean()),
                                 role="descriptive_only" if panel == "challenge8" else "inferential_panel_description"))
    return pd.DataFrame(rows)


def build_outputs(artificial: Path, real: Path) -> dict[str, pd.DataFrame]:
    panels = load_panels(artificial, real)
    return {
        "friedman_current.csv": pd.concat([friedman_rows(p, panels[p]) for p in ("formal21", "real10")], ignore_index=True),
        "pairwise_current.csv": pd.concat([pairwise_rows(p, panels[p]) for p in ("formal21", "real10")], ignore_index=True),
        "descriptive_current.csv": descriptive_rows(panels),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artificial", type=Path, default=ROOT / "results" / "artificial29_six_method_results.csv")
    parser.add_argument("--real", type=Path, default=ROOT / "results" / "real10_six_method_results.csv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs" / "current_statistics")
    args = parser.parse_args()
    protected = {args.artificial.resolve(), args.real.resolve(), (ROOT / "docs" / "benchmark-input-manifest.json").resolve()}
    for filename in ("friedman_current.csv", "pairwise_current.csv", "descriptive_current.csv", "protocol_current.json"):
        if (args.output_dir / filename).resolve() in protected:
            parser.error("An output must not overwrite an input")
    outputs = build_outputs(args.artificial, args.real)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for filename, frame in outputs.items():
        frame.to_csv(args.output_dir / filename, index=False, lineterminator="\n")
    metadata = dict(
        protocol_version=PROTOCOL_VERSION, algorithm_version="v10_0.65", package_version="10.0.65",
        inputs=[dict(role=role, filename=path.name, sha256_canonical_lf=canonical_sha256(path))
                for role, path in (("artificial29_results", args.artificial), ("real10_results", args.real))],
        script_sha256_canonical_lf=canonical_sha256(Path(__file__)),
        input_manifest_sha256_canonical_lf=canonical_sha256(ROOT / "docs" / "benchmark-input-manifest.json"),
        environment=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
        inferential_panels=["formal21", "real10"], descriptive_only=["challenge8"],
        friedman_complete_blocks=dict(formal21=21, real10=9),
        zero_tolerance=EPS, alternative="two-sided", exact_method="conditional sign-permutation dynamic programming",
        holm_family="20 hypotheses (4 metrics x 5 comparators) separately in each inferential panel",
        interval="Hodges-Lehmann Walsh-order 95%; exact no-tie critical rank; conservative with ties; no bootstrap",
        missingness="retain recorded timeout; no quality or count imputation",
        scope="statistics from released scores only; no clustering, label-based metrics, or runtime/RSS rerun",
        outputs={name: dict(rows=len(frame), sha256_canonical_lf=canonical_sha256(args.output_dir / name)) for name, frame in outputs.items()},
    )
    (args.output_dir / "protocol_current.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"Protocol: {PROTOCOL_VERSION}")
    for filename, frame in outputs.items():
        print(f"{filename}: {len(frame)} rows")
    print("formal21/real10 inferential panels; challenge8 descriptive only; no bootstrap or benchmark rerun")


if __name__ == "__main__":
    main()
