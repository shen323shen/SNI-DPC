"""Coverage and missing-value checks for the frozen 39-dataset records."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "verify_release_results.py"
SPEC = importlib.util.spec_from_file_location("verify_release_results", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Cannot load release-result checker")
verify = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(verify)


class ReleaseResultTests(unittest.TestCase):
    def test_frozen_panels(self) -> None:
        self.assertEqual(verify.check_panel(verify.PANELS[0][0], 29), (29, 174, 0))
        self.assertEqual(verify.check_panel(verify.PANELS[1][0], 10), (10, 60, 1))


if __name__ == "__main__":
    unittest.main()
