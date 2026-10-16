# Seven-method statistical records

Protocol: `seven-methods-2026-10-16`. This is a score-level statistical release,
not a new SNI-DPC algorithm revision. Pin the Git commit containing this file;
the checked-in metadata also records the exact input, script, configuration,
and output SHA-256 values.

## Direct access

| Record | Path |
|---|---|
| Artificial inputs: 29 datasets, 203 method records | [`results/artificial29_seven_method_results.csv`](../results/artificial29_seven_method_results.csv) |
| Real inputs: 10 datasets, 70 method records | [`results/real10_seven_method_results.csv`](../results/real10_seven_method_results.csv) |
| All 48 paired tests, effects and intervals | [`results/seven_method_statistics/pairwise_seven.csv`](../results/seven_method_statistics/pairwise_seven.csv) |
| All 740 matched metric differences | [`results/seven_method_statistics/paired_differences_seven.csv`](../results/seven_method_statistics/paired_differences_seven.csv) |
| Eight omnibus tests and SNI-DPC mean ranks | [`results/seven_method_statistics/friedman_seven.csv`](../results/seven_method_statistics/friedman_seven.csv) |
| Available-case means, common-block means and ranks | [`results/seven_method_statistics/descriptive_seven.csv`](../results/seven_method_statistics/descriptive_seven.csv) |
| Input, code and output fingerprints; library versions | [`results/seven_method_statistics/protocol_seven.json`](../results/seven_method_statistics/protocol_seven.json) |
| Evaluation and common configuration record | [`docs/seven-method-evaluation-config.json`](seven-method-evaluation-config.json) |
| Reproducer | [`scripts/reproduce_seven_method_statistics.py`](../scripts/reproduce_seven_method_statistics.py) |

These inputs contain archived scores and actual output counts, not prediction
labels or processed feature matrices. Of 273 method records, 272 are successful;
Gauging-delta on `page_blocks` retains its recorded timeout and empty quality
and count fields. Runtime/RSS fields retain their documented measurement scope.

## Panels and information supplied

The methods are SNI-DPC, Gravity, Torque Clustering, LDCC, Gauging-delta,
DPC-MFP, and WANN-DPC. WANN-DPC receives the integer non-noise reference count,
but not the reference-label vector. The original six receive no reference
count; DPC-MFP uses fixed public defaults rather than automatic neighborhood
parameter selection. Counts in the tables are actual outputs, not the supplied
integers: for example, WANN-DPC produces five groups on `segmentation`, whose
reference count is seven.

All methods retain the registered prepared features and row order. WANN-DPC
bypasses an additional internal Min-Max step. Its reference-label-free
terminal-cycle safeguard is used uniformly; it is triggered on `solar_flare`,
and the other 38 predictions remain pointwise unchanged. The configuration
record documents this call and the common SNI-DPC settings. Data acquisition
and preparation records remain in [the guide](acquisition-and-protocol.md).

- `formal21`: original 21 artificial datasets, an inferential panel.
- `real10`: ten real datasets, a separate inferential panel.
- `challenge8`: eight additional artificial datasets, descriptive only.

The panels are not pooled into a 29- or 39-dataset hypothesis family. Newly
generated configuration-freeze checks are not part of these statistical inputs.

## Omnibus tests and paired comparisons

Higher is better for NMI, ACC, ARI and Purity. Average ranks use mean ranks for
within-dataset ties and rank one for the highest score. All methods use the
same complete blocks: 21 for `formal21`, nine for `real10`. Friedman statistics
include the usual tie correction and use asymptotic chi-square probabilities.
Kendall's W is chi-square divided by `N*(7-1)`. An omnibus result describes a
group difference, not superiority over each comparator.

Paired differences are `SNI-DPC - comparator`, with datasets as the units.
Available matched pairs are retained: 21 for every artificial contrast; ten
for real contrasts except nine against Gauging-delta. The single timeout is
not scored as zero and does not remove `page_blocks` from unrelated pairs.

The paired contrasts are specified in advance and evaluated separately from
the omnibus hypothesis. For each panel:

1. Differences with absolute value at most `1e-12` count as ties and are
   excluded from signed ranks. Tied nonzero magnitudes receive average ranks.
2. The two-sided probability is exact conditional on those observed ranks.
   Twice the ranks are integer weights. Dynamic-programming subset sums
   enumerate the sign-permutation distribution, including both observed tails.
   All-zero differences give probability one and zero positive/negative sums.
3. One Holm family contains **24 hypotheses: four metrics times six
   comparators**. Artificial and real panels are corrected separately.
   The table's `holm_p_joint_24` is the inferential probability, with threshold
   0.05. No within-metric correction substitutes for this family.
4. Rank-biserial effect is `(W_plus-W_minus)/(W_plus+W_minus)`, or zero for
   all-zero differences. Positive values favor SNI-DPC.
5. The location estimate is the Hodges-Lehmann pseudomedian of all Walsh
   averages, including zero paired differences. Its **nominal pointwise 95%**
   interval uses Walsh order statistics and the exact no-tie signed-rank
   critical index. It is not a multiplicity-adjusted simultaneous interval
   or an exact tie-adjusted interval. Small samples or ties affect its
   interpretation; no bootstrap is used.

## Principal real-data findings

Under the panel-wide 24-test Holm correction, positive differences are
significant for NMI against Gravity and DPC-MFP, ACC against DPC-MFP, and ARI
and Purity against Gravity. Each has adjusted probability `0.046875` and
rank-biserial effect `1.0`. For NMI against Gravity, the Hodges-Lehmann shift
is `0.2000065`, with nominal 95% interval `[0.1173532, 0.3229384]` on the
0-to-1 score scale. The other real-data contrasts do not cross the joint
correction threshold; all 24 outcomes, including their directions and
intervals, remain in the complete table.

## Reproduce without clustering

```bash
python scripts/reproduce_seven_method_statistics.py
# Optional independent output directory:
python scripts/reproduce_seven_method_statistics.py --output-dir outputs/seven_check
python -m unittest discover -s tests -p "test_*statistics.py" -v
```

The command validates the frozen input coverage and information conditions,
then writes the four CSVs and one metadata JSON to `outputs/seven_method_statistics/`.
It does not overwrite the checked-in inputs/results, rerun clustering, score
prediction labels, or reproduce runtime/RSS measurements. Compare numeric
results with the checked-in tables; environment fields record the runtime
used for a particular reproduction.

## Preserved earlier entry points

`scripts/reproduce_current_statistics.py` and
`docs/current-statistical-protocol.md` retain the **six-method first-revision**
protocol `current-six-methods-2026-10-15`, with a 20-test Holm family per panel.
Their source and original inputs are unchanged. The seven-method reproducer
reuses its reviewed exact-test and Walsh-interval primitives; metadata pins
that dependency. `scripts/reproduce_statistics.py` remains the separately
labeled earlier v9.5 analysis and is not the seven-method entry point.
