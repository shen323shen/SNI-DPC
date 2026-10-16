# Reproducibility guide

## Environment and software version

The public API is SNI-DPC algorithm `v10_0.65` (Python package version
`10.0.65`). The environment recorded in `requirements-experiments.txt` is
Python 3.10.19, NumPy 2.2.5, SciPy 1.15.3, scikit-learn 1.7.2, pandas
2.3.3 and Matplotlib 3.10.7. The public API needs only NumPy and
scikit-learn; current and historical statistics also require SciPy and pandas.
An environment listing does not by itself establish identical hardware or
execution settings for every archived benchmark record.

## Reproduction levels

1. **Example:** install the package and run `python scripts/run_example.py`.
2. **Equivalence:** run `python -m unittest discover -v`; the public estimator
   is compared pointwise with its checked reference on deterministic inputs.
3. **Record audit:** the seven-method reproducer validates the current 29+10
   coverage, information supplied, status and metric/count fields.
   `python scripts/verify_release_results.py` retains the six-method check.
   `python scripts/verify_checksums.py` checks recorded canonical-LF hashes.
   Neither command reruns experiments or recomputes quality from labels.
4. **Current score-level statistics:** run
   `python scripts/reproduce_seven_method_statistics.py`. It regenerates
   descriptive means/average ranks, Friedman tests, exact two-sided paired
   tests, joint-Holm values, effects and the existing Walsh-order intervals.
   The original 21 and real ten panels remain separate; eight added
   artificial datasets are descriptive only. The real omnibus uses nine
   common complete datasets; paired contrasts use nine or ten as available.
   Each inferential panel has one joint 24-hypothesis Holm family. Full
   tables are checked in under `results/seven_method_statistics/`; the
   [protocol](seven-method-statistical-protocol.md) links every result,
   effect/interval field, configuration record, script and fingerprint.
5. **Full benchmark rerun:** obtain the source distributions and recreate
   the exact processed inputs, acquire compatible comparator implementations,
   and recover the full adapter configurations and stopping rules. The
   [acquisition guide](acquisition-and-protocol.md) records sources, known
   calls, preparation/evaluation rules and unresolved historical details.
   Complete processed inputs, labels and executable six-method runners are
   not bundled; the repository is not a self-contained full rerun package.

For **historical v9.5** only, use
`python scripts/reproduce_statistics.py --output outputs/historical_v9_5`.
This entry point reproduces the old 21-dataset method combination and its
one-sided Wilcoxon/Holm analysis, not the current six-method protocol.
Its calculation rules are preserved rather than silently updated.
The intermediate six-method first-revision entry point is
`scripts/reproduce_current_statistics.py`, protocol
`current-six-methods-2026-10-15`, with 20 joint hypotheses per panel.

## Release scope and input identity

The [release inventory](release-materials.md) distinguishes source code,
score/count tables, runtime/RSS fields, statistical scripts and the example
from unbundled labels, complete processed inputs, raw logs, full execution
configurations and manuscript plotting materials.

The [input manifest](benchmark-input-manifest.json) identifies 39 archived
adapter-input CSVs with source records and raw-byte SHA-256 values.
These are not upstream download hashes or hashes of the standardized
in-memory matrix. A source link or matching shape does not certify exact
reconstruction; unresolved category maps and BNG sampling indices are
explicitly recorded. No private local paths occur in the public manifest.

The SNI-DPC API takes a numeric matrix, not reference labels or `K_true`.
Its internal constants are fixed, so automatic cluster-count selection is not
the same as having no parameters. Five calls are fully automatic;
DPC-MFP uses common fixed defaults without receiving the reference count.
The seventh method, WANN-DPC, is explicitly a known-count control and
receives only the integer non-noise reference count, not the label vector.
Timeout rows have no imputed metrics. Recorded wall-clock time and peak RSS
are study observations, not hardware-independent algorithm properties.
