# Submitted result table

`synthetic21_six_method_results.csv` is the six-method table used for the
submitted 21-dataset artificial benchmark. Its SHA-256 is recorded at the
repository root. The raw column names retain historical provenance; the public
names are mapped as follows:

| Raw prefix | Public method name | True cluster count supplied? |
|---|---|---:|
| `MNG-DPC` | SNI-DPC | No |
| `DPC-oracle-c` | DPC-oracle-c | Yes |
| `WANNDPC` | WANN-DPC | Yes in this evaluated protocol |
| `DSNGCAP` | DSNGCAP | No |
| `Gravity_DSNGCAP` | Gravity+DSNGCAP | No |
| `KNN` | KMeans++ (oracle-c) | Yes |

Run `python scripts/reproduce_statistics.py` to validate the 21 unique datasets
and regenerate the mean NMI/ACC/ARI table and the 15 paired one-sided Wilcoxon
tests with Holm correction. The results describe the evaluated implementations
and protocol; they are not a claim that every comparator universally requires
the same inputs outside this benchmark.

