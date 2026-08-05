from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

from sni_dpc import SNIDPC


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "examples" / "data" / "worked_example_28.csv"
FROZEN = ROOT / "frozen" / "SNI_DPC_v9_5_publication.py"


def load_frozen_module():
    spec = importlib.util.spec_from_file_location("sni_dpc_frozen", FROZEN)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load frozen implementation.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FrozenEquivalenceTests(unittest.TestCase):
    def test_public_api_matches_frozen_partition(self) -> None:
        data = pd.read_csv(DATA).to_numpy(dtype=float)
        frozen_labels, frozen_info = load_frozen_module().mng_dpc(data)
        public = SNIDPC().fit(data)

        np.testing.assert_array_equal(public.labels_, frozen_labels)
        self.assertEqual(public.n_clusters_, frozen_info["num_clusters"])
        self.assertEqual(
            public.natural_neighbor_order_,
            frozen_info["original_info"]["r"],
        )
        np.testing.assert_array_equal(
            public.low_confidence_indices_,
            np.asarray(frozen_info["low_indices"], dtype=int),
        )
