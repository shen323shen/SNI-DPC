from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import unittest

from sni_dpc import SNIDPC


DATA = Path(__file__).resolve().parents[1] / "examples" / "data" / "worked_example_28.csv"
EXPECTED_LABELS = np.asarray([1] * 12 + [0] * 12 + [1, 0, 0, 0], dtype=int)


def load_data() -> np.ndarray:
    return pd.read_csv(DATA).to_numpy(dtype=float)


class WorkedExampleTests(unittest.TestCase):
    def test_worked_example_regression(self) -> None:
        estimator = SNIDPC().fit(load_data())

        np.testing.assert_array_equal(estimator.labels_, EXPECTED_LABELS)
        self.assertEqual(estimator.n_clusters_, 2)
        self.assertEqual(estimator.natural_neighbor_order_, 6)
        np.testing.assert_array_equal(
            estimator.low_confidence_indices_,
            np.asarray([9, 10, 22, 24, 25, 26, 27]),
        )

    def test_invalid_input_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "two-dimensional"):
            SNIDPC().fit([1, 2, 3])

        with self.assertRaisesRegex(ValueError, "NaN"):
            SNIDPC().fit(np.asarray([[0.0, np.nan], [1.0, 1.0]]))
