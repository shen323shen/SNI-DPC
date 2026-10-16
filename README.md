# SNI-DPC

**SNI-DPC: A Screened Natural-Neighbor Interaction Field for Automatic Density
Peaks Clustering**

This repository contains the SNI-DPC `v10_0.65` implementation and audited
dataset-level results for 29 artificial and 10 real-world datasets. The
method uses an extended natural-neighbor relation and a screened structural
field force to assemble clusters and determine their number without supplying
the true cluster count. The fixed structural-field density weight is 0.65;
the other three density views share the remaining 0.35. These are internal
constants, not a claim that the method has no parameters.

The algorithm revision is named `v10_0.65`; the installable Python package
uses the PEP 440 version `10.0.65`.

## Quick start

The public API requires Python 3.10--3.12, NumPy, and scikit-learn.

```bash
python -m pip install -e ".[test,experiments]"
python scripts/run_example.py
python -m unittest discover -v
```

The example is a fixed 28-point input. It should report `2` clusters, natural-
neighbor order `6`, and `7` low-confidence points. The true cluster count is
not passed to the algorithm.

```python
import numpy as np
from sni_dpc import SNIDPC

X = np.loadtxt("examples/data/worked_example_28.csv", delimiter=",", skiprows=1)
model = SNIDPC().fit(X)
labels = model.labels_
print(model.n_clusters_, model.natural_neighbor_order_)
```

## Repository map

| Path | Purpose |
|---|---|
| `src/sni_dpc/` | Public API and self-contained v10_0.65 computational core |
| `frozen/` | Checked v10_0.65 reference and preserved v9.5 historical snapshot |
| `examples/` | Small deterministic input for a complete end-to-end run |
| `results/` | Current seven-method inputs/statistics, preserved six-method records and earlier historical tables |
| `scripts/` | Example, integrity checks, current statistics and a separately labeled historical entry point |
| `tests/` | Public-API equivalence and table-integrity tests |
| `docs/` | Implementation map, data provenance, and reproducibility notes |

## Current result records

The seven-method comparison adds **WANN-DPC as a known-count control** to
SNI-DPC, Gravity, Torque Clustering, LDCC, Gauging-delta and DPC-MFP.
`results/artificial29_seven_method_results.csv` and
`results/real10_seven_method_results.csv` contain 273 records: 272 successes
and the retained `page_blocks`/Gauging-delta timeout. Only WANN-DPC receives
the integer non-noise reference count. It uses the common prepared features
without an additional internal Min-Max transformation and the documented
reference-label-free cycle safeguard.

The complete corrected paired tests, effects, nominal intervals, omnibus
results, configuration record and code fingerprints are directly available
under [`results/seven_method_statistics/`](results/seven_method_statistics/).
See the [seven-method protocol](docs/seven-method-statistical-protocol.md)
and [evaluation configuration](docs/seven-method-evaluation-config.json).

```bash
python scripts/reproduce_seven_method_statistics.py
```

Protocol `seven-methods-2026-10-16` keeps the original 21 artificial and ten
real datasets as separate inferential panels, with 24 panel-wide Holm tests
each; the eight additional artificial datasets remain descriptive. This is
statistics from archived scores, not a clustering or resource rerun.

### Preserved six-method first-revision records

`results/artificial29_six_method_results.csv` and
`results/real10_six_method_results.csv` record one row per dataset and method
for SNI-DPC, Gravity, Torque Clustering, LDCC, Gauging-delta and DPC-MFP.
Run `python scripts/verify_release_results.py` to check coverage, run status
and metric/count fields. The `page_blocks`/Gauging-delta run timed out; no
quality metrics are imputed for it.

Run `python scripts/reproduce_current_statistics.py` to reproduce six-method statistics
from the published scores: the original 21 artificial datasets and ten real
datasets have separate inferential panels; the eight artificial extensions
remain descriptive. Outputs go to `outputs/current_statistics/`. See the
[preserved six-method protocol](docs/current-statistical-protocol.md) for exact
two-sided signed-rank tests, panel-specific Holm families, effects and the
existing non-bootstrap intervals. This command does not rerun clustering or
recompute metrics from prediction labels.

These are archived experiment records, not a benchmark rerun by this
repository. Their hashes and reproduction boundaries are described in
`checksums.sha256` and `docs/reproducibility.md`. The former 21-dataset table
and `scripts/reproduce_statistics.py` remain available as **v9.5 history**;
they describe neither the preserved six-method nor the current seven-method comparison.

## Reproducibility boundaries

The public API is regression-tested against the checked v10_0.65 reference on
deterministic inputs. Full benchmark reruns require the executed dataset
versions, preprocessing and experiment configuration. Only SNI-DPC source is
distributed here; third-party comparator code and input data without verified
redistribution terms are not included.

Prediction labels, the complete processed-input archive, raw monitoring logs,
the complete executable comparator runners and the full figure-generation
archive are not released. The [material inventory](docs/release-materials.md)
distinguishes them from the source code, result fields, example and scripts
that are included. The [39-input identity manifest](docs/benchmark-input-manifest.json)
provides source records and archived CSV SHA-256 values without exposing
private filesystem paths or distributing inputs.

See the [acquisition and preparation guide](docs/acquisition-and-protocol.md),
[data provenance](docs/data-provenance.md), and
[reproducibility guide](docs/reproducibility.md) for access steps and the
remaining requirements for a full benchmark rerun.

## Citation

Use the metadata in `CITATION.cff`. The corresponding author is Lijun Yang
(`ylijun@swun.edu.cn`).

## License

The SNI-DPC source in this repository is released under the BSD 3-Clause
License. Dependencies retain their own licenses.

