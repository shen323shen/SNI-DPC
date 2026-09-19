# Reproducibility guide

## Environment and software version

The public API is SNI-DPC algorithm `v10_0.65` (Python package version
`10.0.65`). The environment recorded in `requirements-experiments.txt` is
Python 3.10.19, NumPy 2.2.5, SciPy 1.15.3, scikit-learn 1.7.2, pandas
2.3.3 and Matplotlib 3.10.7. The public API needs only NumPy and
scikit-learn; historical summary reproduction also requires SciPy and pandas.
An environment listing does not by itself establish identical hardware or
execution settings for every archived benchmark record.

## Reproduction levels

1. **Example:** install the package and run `python scripts/run_example.py`.
2. **Equivalence:** run `python -m unittest discover -v`; the public v10_0.65
   estimator is compared pointwise with its checked reference on several
   deterministic inputs.
3. **Record audit:** run `python scripts/verify_release_results.py`. It checks
   the 29+10 dataset coverage, six methods, status and metric/count fields;
   run `python scripts/verify_checksums.py` to check the recorded canonical-LF
   content hashes on any platform. It does not
   rerun experiments or recompute quality metrics from prediction labels.
4. **Full rerun:** acquire the executed input versions and experiment
   manifests, reproduce the recorded preprocessing, comparator implementations
   and stopping rules, then run the complete comparison. Input datasets and
   third-party baseline sources are not bundled here.

For the **historical v9.5** 21-dataset comparison only, run
`python scripts/reproduce_statistics.py` to regenerate its mean metrics and
paired Wilcoxon/Holm summaries. Those statistics are not tests of the current
39-dataset six-method results.

The SNI-DPC API takes a numeric matrix, not reference labels or `K_true`.
Its internal constants are fixed, so automatic cluster-count selection is not
the same as having no parameters. Timeout rows have no imputed metrics.
Wall-clock time and peak RSS reflect recorded runs rather than
hardware-independent algorithm properties.
