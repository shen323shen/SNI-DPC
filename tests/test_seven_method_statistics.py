"""Seven-method protocol and input guards; the six-method regression stays intact."""
import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location("seven_statistics", ROOT / "scripts/reproduce_seven_method_statistics.py")
seven = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(seven)
ART = ROOT / "results/artificial29_seven_method_results.csv"
REAL = ROOT / "results/real10_seven_method_results.csv"


class SevenMethodTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.panels = seven.load_panels(ART, REAL)
        cls.outputs = seven.build_outputs(ART, REAL)

    def test_panel_sizes_and_reference_count_access(self):
        self.assertEqual({k: len(v) for k, v in self.panels.items()}, {"formal21": 147, "challenge8": 56, "real10": 70})
        for frame in self.panels.values():
            self.assertTrue(frame.reference_count_supplied.equals(frame.method.eq("WANN-DPC")))

    def test_complete_omnibus_blocks(self):
        frame = self.outputs["friedman_seven.csv"]
        self.assertEqual(len(frame), 8)
        self.assertEqual(set(frame.methods), {7})
        self.assertEqual(set(frame[frame.panel == "formal21"].complete_datasets), {21})
        self.assertEqual(set(frame[frame.panel == "real10"].complete_datasets), {9})
        self.assertTrue(np.allclose(frame.kendall_w, frame.friedman_chi_square / (frame.complete_datasets * 6)))

    def test_panel_wide_holm24_and_pair_counts(self):
        frame = self.outputs["pairwise_seven.csv"]
        self.assertEqual(frame.groupby("panel").size().to_dict(), {"formal21": 24, "real10": 24})
        for _, part in frame.groupby("panel"):
            np.testing.assert_allclose(part.holm_p_joint_24, seven.holm_adjust(part.wilcoxon_p_two_sided_exact_signflip), rtol=0, atol=0)
        for row in frame.itertuples():
            expected = 21 if row.panel == "formal21" else 9 if row.comparator == "Gauging-delta" else 10
            self.assertEqual(row.paired_datasets, expected)
            self.assertEqual(row.wins + row.ties + row.losses, expected)

    def test_all_differences_preserved_and_challenges_descriptive(self):
        self.assertEqual(len(self.outputs["paired_differences_seven.csv"]), 740)
        self.assertEqual(len(self.outputs["descriptive_seven.csv"]), 84)
        self.assertNotIn("challenge8", set(self.outputs["pairwise_seven.csv"].panel))

    def test_actual_wann_count_is_not_replaced_by_given_count(self):
        row = self.panels["real10"].query("dataset=='segmentation' and method=='WANN-DPC'").iloc[0]
        self.assertEqual((row.true_c, row.c_out), (7, 5))

    def test_duplicate_rejected(self):
        frame = pd.read_csv(ART)
        with patch.object(seven.pd, "read_csv", return_value=pd.concat([frame, frame.iloc[:1]])):
            with self.assertRaises(ValueError):
                seven.load_panels(ART, REAL)

    def test_timeout_cannot_be_zero_score(self):
        frame = pd.read_csv(REAL)
        frame.loc[frame.status == "timeout", "NMI"] = 0
        with patch.object(seven.pd, "read_csv", side_effect=[pd.read_csv(ART), frame]):
            with self.assertRaises(ValueError):
                seven.load_panels(ART, REAL)

    def test_no_count_access_for_other_methods(self):
        frame = pd.read_csv(ART)
        frame.loc[frame.method == "SNI-DPC", "reference_count_supplied"] = True
        with patch.object(seven.pd, "read_csv", return_value=frame):
            with self.assertRaises(ValueError):
                seven.load_panels(ART, REAL)


if __name__ == "__main__":
    unittest.main()
