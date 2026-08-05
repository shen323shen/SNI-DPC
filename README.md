# SNI-DPC

**SNI-DPC: A Screened Natural-Neighbor Interaction Field for Automatic Density
Peaks Clustering**

This repository contains the public SNI-DPC implementation, a deterministic
worked example, and the result table used for the submitted Knowledge-Based
Systems manuscript. SNI-DPC constructs a screened interaction field on a
natural-neighbor graph and discovers both the natural-neighbor order and the
output cluster count from the data.

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
| `src/sni_dpc/` | Clean public API and the extracted computational core |
| `frozen/` | Immutable publication implementation used for the submitted results |
| `examples/` | Small deterministic input for a complete end-to-end run |
| `results/` | Formal 21-dataset six-method table and regenerated summaries |
| `scripts/` | Core extraction, example execution, and statistical reproduction |
| `tests/` | Standard-library regression tests against the frozen implementation |
| `docs/` | Implementation map, data provenance, and reproducibility notes |

## Reproducing the paper evidence

The formal table is `results/synthetic21_six_method_results.csv`. It contains
21 artificial datasets and six evaluated methods. Historical raw column names
are retained for traceability; the public names and true-cluster-count policy
are documented in `results/README.md`.

Regenerate the aggregate metrics and all 15 paired one-sided Wilcoxon tests with
Holm correction:

```bash
python scripts/reproduce_statistics.py
```

The script writes `results/mean_metrics_21.csv` and
`results/wilcoxon_holm_21.csv`. These are summaries of the frozen result table;
they do not silently rerun or alter the submitted experiments.

## Reproducibility boundaries

The frozen implementation and formal result table are checksum protected. The
clean API is tested against the frozen implementation on the deterministic
worked example. Full benchmark reruns require the original dataset files and
the recorded execution environment; the repository does not redistribute data
whose exact historical provenance or license is unresolved.

The repository contains SNI-DPC code only. DSNGCAP and Gravity+DSNGCAP are
reported as evaluated comparator modules in the paper and result table, but
their source code is not redistributed here because their publication and
group-ownership terms must be handled separately. The same applies to any
third-party implementation not covered by this repository license.

See `docs/data-provenance.md` and `docs/reproducibility.md` for the exact scope.

## Citation

Use the metadata in `CITATION.cff`. The corresponding author is Lijun Yang
(`ylijun@swun.edu.cn`).

## License

The SNI-DPC source in this repository is released under the BSD 3-Clause
License. Dependencies retain their own licenses.

