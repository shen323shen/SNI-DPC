# Current six-method statistical protocol

Protocol version: `current-six-methods-2026-10-15`.
Algorithm: SNI-DPC `v10_0.65`; Python package: `10.0.65`.

## Inputs and panels

The script reads the two frozen public result CSVs, not prediction labels.
The six method columns are SNI-DPC, Gravity, Torque Clustering, LDCC,
Gauging-delta, and DPC-MFP. No current call receives the reference count;
DPC-MFP is a fixed-default call, distinct from the other five automatic calls.

- `formal21`: the original 21 artificial datasets, the primary inferential panel.
- `challenge8`: the eight added artificial datasets, descriptive only.
- `real10`: ten real datasets, a separate inferential panel. The public source
  CSV retains the old collection name `real_top30`; the reader explicitly
  maps this name to `real10` without modifying the source file.

The script validates dataset/method uniqueness, membership against the public
identity manifest, successful metric ranges and the one recorded
`page_blocks`/Gauging-delta timeout. Missing quality/count fields are never
replaced by zero. The panels are not pooled to create a new 39-dataset test.

## Descriptive summaries and omnibus tests

NMI, ACC, ARI and Purity are oriented so that higher is better. Descriptive
means use each method's available successful datasets and state that count.
Average ranks use the same complete six-method blocks for all six methods;
within a dataset, tied scores receive their mean rank and rank 1 is best.

Friedman tests use 21 common complete blocks for `formal21`, and nine for
`real10`. The latter excludes `page_blocks` only from the common-block
omnibus calculation, because Gauging-delta did not complete. Kendall's
concordance is `chi_square_F / (N * (6 - 1))`. An omnibus difference does not
by itself establish SNI-DPC superiority.

## Paired tests, multiplicity and effects

The paired unit is a dataset. Differences are `SNI-DPC - comparator`.

1. Differences with `abs(d) <= 1e-12` remain ties in win/tie/loss counts and
   are excluded from signed ranks. Nonzero tied magnitudes receive average ranks.
2. The probability is **two-sided and exact conditional on the observed
   ranks**. Twice each rank is an integer; dynamic-programming subset sums
   enumerate the sign-permutation distribution. Both inclusive tails are
   calculated and `p = min(1, 2 * min(lower, upper))`. No normal approximation
   or SciPy `method="auto"` selection is used. All-zero differences give
   `p=1`, `W_plus=W_minus=0`.
3. One joint Holm correction covers 20 hypotheses (four metrics times five
   comparators) **separately in each panel**. These are not a 40-test family.
   The five-test within-metric Holm column is diagnostic only.
4. Rank-biserial effect is `(W_plus-W_minus)/(W_plus+W_minus)`, or zero for
   an all-zero comparison. The location estimate is the Hodges--Lehmann
   pseudomedian of all Walsh averages, including zero paired differences.
5. The existing 95% interval uses Walsh order statistics with the exact
   no-tie signed-rank critical index. With ties or small panels it can be
   conservative; it is not a bootstrap interval and is not described as an
   exact tie-adjusted confidence interval. No new bootstrap is performed.

Real-data paired comparisons use ten available pairs for Gravity, Torque
Clustering, LDCC and DPC-MFP, and nine for Gauging-delta. The extension panel
does not enter either inferential family. Statistical significance requires
the joint-Holm value below 0.05; numerical advantages are reported separately.

## Command, outputs and provenance

```bash
python scripts/reproduce_current_statistics.py
# Optional: keep outputs in a separate directory.
python scripts/reproduce_current_statistics.py --output-dir outputs/current_check
```

- `friedman_current.csv`: eight panel/metric omnibus records.
- `pairwise_current.csv`: 40 records, with exact probabilities, counts,
  effects, existing intervals and Holm corrections.
- `descriptive_current.csv`: 72 panel/metric/method descriptions.
- `protocol_current.json`: protocol and algorithm versions, input/script
  canonical-LF SHA-256 values, local library versions, output checksums and scope.

These outputs reproduce statistics from the archived scores. They do not
recompute the scores, runtime or memory. The frozen input tables remain unchanged.

## Historical entry point

`scripts/reproduce_statistics.py` belongs to **historical v9.5** and the older
21-artificial-dataset comparison, including oracle-count baselines. It uses
one-sided Wilcoxon tests with SciPy's automatic calculation selection.
It is preserved for historical reproduction and is incompatible with the
current two-sided protocol. Use a distinct output directory:

```bash
python scripts/reproduce_statistics.py --output outputs/historical_v9_5
```
