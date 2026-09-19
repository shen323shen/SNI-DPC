# Benchmark result records

## Current v10_0.65 comparison

`artificial29_six_method_results.csv` and `real10_six_method_results.csv`
contain 29 artificial and 10 real datasets, respectively. Each dataset has
one record for each of SNI-DPC, Gravity, Torque Clustering, LDCC,
Gauging-delta and DPC-MFP. All six methods belong to the current automatic
cluster-count comparison; none is presented as an oracle-count baseline.
The files retain the source tables' original column names, ordering and
numerical precision. `true_c` is the reference cluster count for evaluation;
`c_out` is the predicted cluster count. The record for
`page_blocks`/Gauging-delta has `status=timeout`, blank quality/count values
and elapsed runtime of approximately 8.2 hours. Missing values are not zeros.

These files were transferred without modifying records or numeric precision
from the reviewed six-method tables
`R3_05_ARTIFICIAL29_SIX_METHOD_RESULTS_planB.csv` and
`R2_05_REAL_NATIVE_SIX_METHOD_RESULTS_planB.csv`. Run
`python scripts/verify_release_results.py` for structural validation and
run `python scripts/verify_checksums.py` to check canonical-LF content hashes
on any platform. The validator does not rerun experiments. Preprocessing,
metric conventions and timeout policies
must be read with the executed experiment manifests, not inferred from the
CSV alone.

## Historical v9.5 submission (not the current benchmark)

`synthetic21_six_method_results.csv`, `mean_metrics_21.csv` and
`wilcoxon_holm_21.csv` preserve the earlier 21-artificial-dataset comparison.
The historical raw column mapping was:

| Raw prefix | Public method name | True cluster count supplied? |
|---|---|---:|
| `MNG-DPC` | SNI-DPC v9.5 | No |
| `DPC-oracle-c` | DPC-oracle-c | Yes |
| `WANNDPC` | WANN-DPC | Yes in this evaluated protocol |
| `DSNGCAP` | DSNGCAP | No |
| `Gravity_DSNGCAP` | Gravity+DSNGCAP | No |
| `KNN` | KMeans++ (oracle-c) | Yes |

`python scripts/reproduce_statistics.py` regenerates only the historical mean
NMI/ACC/ARI table and 15 paired Wilcoxon/Holm tests. These historical
oracle-count baselines must not be mixed with the current six-method automatic
comparison.
