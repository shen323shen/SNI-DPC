# Released materials and access scope

Inventory date: 16 October 2026. The algorithm is SNI-DPC `v10_0.65`, package
`10.0.65`, unchanged by this statistical release. Current statistics use
`seven-methods-2026-10-16`; the earlier six-method protocol is preserved.
Use a Git commit to identify the repository snapshot; these are distinct
version identifiers, not a new algorithm release.

| Material | In this repository | What the release supports |
|---|---|---|
| SNI-DPC source and fixed internal settings | Yes: `src/sni_dpc/`, `frozen/` | Run the method and compare the public estimator with the checked reference |
| Current seven-method score/count results | Yes: `results/artificial29_seven_method_results.csv`, `results/real10_seven_method_results.csv` | Inspect 273 dataset-method records, including 272 successes and one timeout; WANN-DPC is a known-count control |
| Runtime and peak-RSS observations | Fields in the result tables, with resource-scope notes | Inspect archived six-method observations and independent WANN real-data timings; unmeasured fields remain empty |
| Current statistical reproduction | Yes: `scripts/reproduce_seven_method_statistics.py` and `docs/seven-method-statistical-protocol.md` | Recompute panel summaries, Friedman tests, 24-test-per-panel Holm adjustments, paired effects and nominal Walsh intervals |
| Complete current statistical tables | Yes: `results/seven_method_statistics/` | Directly access 48 paired tests, 740 metric differences, 8 omnibus results, 84 descriptions and checksummed metadata |
| Evaluation/configuration record | Yes: `docs/seven-method-evaluation-config.json` | Inspect information supplied, common SNI-DPC controls, prepared-input rules, WANN adaptation and executed-source fingerprints; not a complete third-party runner |
| Preserved six-method first-revision records | Original two six-method CSVs, `scripts/reproduce_current_statistics.py` and its protocol | Reproduce the unchanged 20-test-per-panel analysis separately from the current seven-method release |
| Prediction labels for the 39-dataset comparison | No | Scores cannot be recomputed from predictions using this release alone |
| Complete processed benchmark inputs | No | Sources and archived input identities are documented, not a bundled input archive |
| Dataset identity manifest | Yes: `docs/benchmark-input-manifest.json` | Inspect 39 names, panels, sizes, dimensions, reference counts, source records and archived CSV SHA-256 values |
| Acquisition and preprocessing instructions | Yes: `docs/data-provenance.md`, `docs/acquisition-and-protocol.md` | Locate upstream distributions and follow the recorded preparation/evaluation rules; unresolved historical encoding/sample details remain identified |
| Complete executable comparator runners | No | Documented method calls and common settings are provided, but the complete executable benchmark workflow is not bundled |
| Third-party method implementations | No | Obtain separately from the cited upstream repositories or publication records; check their licenses and runtime requirements |
| Complete figure-generation materials | No | The illustrative example is runnable; the complete manuscript plotting archive is not included |
| Historical v9.5 scores and statistical script | Yes, explicitly historical | Reproduce only the older 21-dataset method combination, not the current seven-method analysis |
| Small end-to-end example | Yes: `examples/data/worked_example_28.csv` and `scripts/run_example.py` | Run the 28-point example without any benchmark dataset download |

## What can be reproduced directly

Install the package, run the example and regression tests, validate the
published records, and recompute current statistics:

```bash
python -m pip install ".[test,experiments]"
python scripts/run_example.py
python -m unittest discover -v
python scripts/verify_checksums.py
python scripts/verify_release_results.py
python scripts/reproduce_seven_method_statistics.py
# Preserved six-method first-revision statistics:
python scripts/reproduce_current_statistics.py
```

Generated files go to the ignored `outputs/seven_method_statistics/` or
`outputs/current_statistics/`, not over the checked-in result tables. No labels,
input data, or clustering runs are needed for this score-level reproduction.

## Full benchmark reruns

A full rerun additionally requires the source data, the exact processed
representations, third-party implementations and the complete execution
configuration. Start with the per-dataset source records and method-specific
access instructions in [Acquisition and protocol](acquisition-and-protocol.md).
The identity manifest identifies archived numeric CSV files; it does not
distribute those files or certify a byte-identical reconstruction from each
upstream download. Historical sampling indices, category mappings and full
adapter configurations that are not published remain requirements for an
exact rerun.

The published runtime and RSS values are archived observations. Recomputing
statistics from those tables does not reproduce resource measurements, the
original hardware workload, or metrics calculated from prediction labels.
