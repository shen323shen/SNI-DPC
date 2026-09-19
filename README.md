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
python -m pip install -e .
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
| `results/` | 29-artificial/10-real six-method records and historical v9.5 tables |
| `scripts/` | Example, core-generation checks, result validation and historical summaries |
| `tests/` | Public-API equivalence and table-integrity tests |
| `docs/` | Implementation map, data provenance, and reproducibility notes |

## Current result records

`results/artificial29_six_method_results.csv` and
`results/real10_six_method_results.csv` record one row per dataset and method
for SNI-DPC, Gravity, Torque Clustering, LDCC, Gauging-delta and DPC-MFP.
Run `python scripts/verify_release_results.py` to check coverage, run status
and metric/count fields. The `page_blocks`/Gauging-delta run timed out; no
quality metrics are imputed for it.

These are archived experiment records, not a benchmark rerun by this
repository. Their hashes and reproduction boundaries are described in
`checksums.sha256` and `docs/reproducibility.md`. The former 21-dataset table
and `scripts/reproduce_statistics.py` remain available as **v9.5 history**;
they do not describe the current six-method comparison.

## Reproducibility boundaries

The public API is regression-tested against the checked v10_0.65 reference on
deterministic inputs. Full benchmark reruns require the executed dataset
versions, preprocessing and experiment configuration. Only SNI-DPC source is
distributed here; third-party comparator code and input data without verified
redistribution terms are not included.

See `docs/data-provenance.md` and `docs/reproducibility.md` for the exact scope.

## Citation

Use the metadata in `CITATION.cff`. The corresponding author is Lijun Yang
(`ylijun@swun.edu.cn`).

## License

The SNI-DPC source in this repository is released under the BSD 3-Clause
License. Dependencies retain their own licenses.

