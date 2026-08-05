# Reproducibility guide

## Environment

The submitted experiments used Python 3.10.19 with NumPy 2.2.5, SciPy 1.15.3,
scikit-learn 1.7.2, pandas 2.3.3, and Matplotlib 3.10.7. The exact package
records are listed in `requirements-experiments.txt`. The public clustering API
needs only NumPy and scikit-learn.

## Reproduction levels

1. **Smoke run:** install the package and run `python scripts/run_example.py`.
2. **Regression run:** execute `python -m unittest discover -v`. This compares
   the public API with the frozen publication implementation on the deterministic
   28-point input.
3. **Table audit:** run `python scripts/reproduce_statistics.py`. This checks
   the formal table, recomputes six-method means, and recomputes 15 paired
   one-sided Wilcoxon tests with Holm correction.
4. **Full experiment rerun:** obtain the source data listed in the manuscript
   and Supplementary Material, record the exact files and preprocessing in a
   manifest, and use the frozen implementation. This level is intentionally not
   part of the default CI workflow because some inherited dataset versions are
   not recoverable from current public package metadata.

## What is guaranteed

- SNI-DPC does not receive the true number of clusters in the public API.
- The worked example output is deterministic and regression-tested.
- The formal CSV and frozen implementation have immutable SHA-256 values.
- The statistical script does not modify the formal table.

## What is not claimed

Runtime values describe the evaluated implementation and workstation. The
eight real-world datasets are a bounded external-validity study, not a claim of
universal real-data superiority. A completed single-cluster output with zero
NMI/ARI is a valid result unless the run status proves failure.

