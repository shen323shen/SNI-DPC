"""Current protocol checks; intentionally independent of historical v9.5 tests."""

import importlib.util
import itertools
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("current_statistics", ROOT / "scripts" / "reproduce_current_statistics.py")
current = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(current)
ARTIFICIAL = ROOT / "results" / "artificial29_six_method_results.csv"
REAL = ROOT / "results" / "real10_six_method_results.csv"


class ExactTests(unittest.TestCase):
    def test_all_zero(self):
        self.assertEqual(current.exact_signed_rank_pvalue(np.array([0, 1e-13, -1e-12])), (1.0, 0.0, 0.0))

    def test_ties_match_brute_force(self):
        for values in ([1, -1, 2, 3, 0], [1, 2, 3, 4], [-1, -2, -3, -4], [1, -1, 1, -1]):
            values = np.asarray(values, dtype=float)
            nz = values[np.abs(values) > current.EPS]
            ranks = stats.rankdata(np.abs(nz), method="average")
            observed = ranks[nz > 0].sum()
            sums = [ranks[np.asarray(mask, dtype=bool)].sum() for mask in itertools.product((0, 1), repeat=len(nz))]
            p = min(1, 2 * min(np.mean(np.asarray(sums) <= observed), np.mean(np.asarray(sums) >= observed)))
            self.assertEqual(current.exact_signed_rank_pvalue(values)[0], p)

    def test_sign_reversal(self):
        a = current.exact_signed_rank_pvalue(np.array([1, -2, 3, 3, -4]))
        b = current.exact_signed_rank_pvalue(np.array([-1, 2, -3, -3, 4]))
        self.assertEqual(a, (b[0], b[2], b[1]))

    def test_nonfinite_rejected(self):
        with self.assertRaises(ValueError):
            current.exact_signed_rank_pvalue(np.array([float("nan")]))

    def test_holm(self):
        np.testing.assert_allclose(current.holm_adjust([0.01, 0.04, 0.03, 0.2]), [0.04, 0.09, 0.09, 0.2])

    def test_holm_invalid(self):
        with self.assertRaises(ValueError):
            current.holm_adjust([float("nan")])

    def test_walsh_all_zero(self):
        self.assertEqual(current.hodges_lehmann_interval(np.zeros(10))[:3], (0.0, 0.0, 0.0))

    def test_walsh_tiny_panel(self):
        hl, low, high, critical = current.hodges_lehmann_interval(np.array([1.0]))
        self.assertEqual((hl, critical), (1.0, -1))
        self.assertTrue(np.isnan(low) and np.isnan(high))


class CurrentRecordsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.panels = current.load_panels(ARTIFICIAL, REAL)
        cls.outputs = current.build_outputs(ARTIFICIAL, REAL)

    def test_panel_coverage(self):
        self.assertEqual({p: len(f) for p, f in self.panels.items()}, dict(formal21=126, challenge8=48, real10=60))

    def test_complete_blocks(self):
        friedman = self.outputs["friedman_current.csv"]
        self.assertEqual(set(friedman[friedman.panel == "formal21"].complete_datasets), {21})
        self.assertEqual(set(friedman[friedman.panel == "real10"].complete_datasets), {9})
        self.assertNotIn("challenge8", set(friedman.panel))

    def test_reference_friedman_values(self):
        expected = [44.24722662440571, 46.0179153094462, 42.66558966074309, 2.0506912442396463,
                    26.942675159235645, 9.426751592356682, 24.600638977635768, 30.657051282051256]
        # Reviewed fixed statistics, in NMI/ACC/ARI/Purity order for formal21 then real10.
        np.testing.assert_allclose(self.outputs["friedman_current.csv"].friedman_chi_square, expected, rtol=0, atol=1e-12)

    def test_pair_counts_and_holm_families(self):
        tests = self.outputs["pairwise_current.csv"]
        self.assertEqual(tests.groupby("panel").size().to_dict(), dict(formal21=20, real10=20))
        for _, row in tests.iterrows():
            self.assertEqual(row.paired_datasets, 9 if row.panel == "real10" and row.comparator == "Gauging-delta" else 10 if row.panel == "real10" else 21)
            self.assertEqual(row.wins + row.ties + row.losses, row.paired_datasets)
        for _, group in tests.groupby("panel"):
            np.testing.assert_allclose(group.holm_p_joint_20, current.holm_adjust(group.wilcoxon_p_two_sided_exact_signflip))

    def test_manifest_private_paths_absent(self):
        text = (ROOT / "docs" / "benchmark-input-manifest.json").read_text(encoding="utf-8")
        manifest = json.loads(text)
        self.assertEqual(len(manifest["datasets"]), 39)
        self.assertNotIn("dataset_path", text)
        self.assertNotIn("Users\\", text)

    def test_duplicate_rejected(self):
        frame = pd.read_csv(ARTIFICIAL)
        with patch.object(current.pd, "read_csv", return_value=pd.concat([frame, frame.iloc[:1]])):
            with self.assertRaises(ValueError):
                current.load_panels(ARTIFICIAL, REAL)

    def test_timeout_zero_imputation_rejected(self):
        frame = pd.read_csv(REAL)
        frame.loc[frame.status == "timeout", "NMI"] = 0
        artificial = pd.read_csv(ARTIFICIAL)
        with patch.object(current.pd, "read_csv", side_effect=[artificial, frame]):
            with self.assertRaises(ValueError):
                current.load_panels(ARTIFICIAL, REAL)

    def test_history_not_accepted(self):
        with self.assertRaises(ValueError):
            current.load_panels(ROOT / "results" / "synthetic21_six_method_results.csv", REAL)


if __name__ == "__main__":
    unittest.main()
