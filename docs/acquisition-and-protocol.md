# Dataset acquisition and recorded benchmark protocol

This guide separates upstream acquisition, preparation rules and the archived
numeric inputs. It documents the current 29-artificial/10-real comparison;
it is not a complete executable comparator benchmark runner. The current
seven-method statistical release is specified in
[the seven-method protocol](seven-method-statistical-protocol.md).

## 1. Locate the correct source distribution

Each entry in [benchmark-input-manifest.json](benchmark-input-manifest.json)
contains the public result-table name, panel, sample count, feature count,
non-noise reference count, source link and archived adapter-input CSV SHA-256.
Follow the source details in [data-provenance.md](data-provenance.md).

- **OpenML:** open the `/d/<data_id>` record in the manifest, verify its
  version, target field and (where recorded) file ID, and use its download
  link or API. The metadata endpoint is
  `https://www.openml.org/api/v1/json/data/<data_id>`; data are downloadable
  from `https://www.openml.org/data/v1/download/<file_id>` when a file ID is
  recorded. Do not substitute a same-name upload with another ID/version.
- **UCI:** follow the persistent DOI for banknote or Image Segmentation and
  download the record's source files. The evaluated segmentation table
  appends `segmentation.test` after `segmentation.data`. UCI source downloads
  and the processed numeric table have different byte hashes.
- **Artificial benchmarks:** follow the linked source family/artifact. For
  immutable ARFF/MAT links, use that commit's raw file, preserve row order,
  and parse features and reference labels separately. A source family or
  generator command is not proof that a new draw reproduces an inherited file.

No prediction labels, full processed-input archive or full frozen adapter
configuration is bundled. Original BNG sample indices and historical category
maps are not supplied by this release. Matching source ID, `n`, `d` or label
count alone is insufficient to recover exact rows and encoding; those cases
remain qualified rather than represented as exact reconstructions.

## 2. Preparation and label isolation

The archived numeric format has `d` feature columns followed by one reference
label column. Separate the label column **before** fitting. The original six
calls receive neither reference labels nor the reference count. The added
WANN-DPC control receives only the integer number of non-noise reference
classes, not the label vector. Preserve sample order and keep evaluation
labels aligned outside each method.

- The original 21 artificial inputs retain their numerical coordinates;
  `hypercube` uses all three dimensions.
- The added artificial inputs retain the registered numerical representation.
  BNG inputs are fixed 600-row processed samples, not the full OpenML records.
  The evaluated file named `2d-10c` has **nine nonempty reference labels**,
  despite the ten-class declaration in its inherited source header.
- The ten real cases use fixed categorical encoding and missing-value
  treatment, followed by per-feature z-score standardization. Constant
  columns/non-finite standardized results become zero. The full historical
  category order and sample-level transformation manifests are not released;
  applying these broad rules does not certify byte-identical inputs.
- Original six-method internal normalization remains part of each implementation.
  WANN-DPC uses the same prepared matrix without an additional internal
  Min-Max step; its registered no-reference-label cycle safeguard is retained.
  PCA belongs to the separately controlled dimension experiments and is not
  applied silently to the native-dimension ten-real benchmark.

The manifest SHA-256 hashes identify the archived adapter-input CSV bytes,
before any adapter-side normalization. They do not hash the upstream ARFF/MAT
or the in-memory standardized matrix. CSV serialization, headers, column order
and row order affect this hash even when numeric values are equivalent.

## 3. Obtain comparators and retain the documented calls

Only SNI-DPC code is distributed here. Third-party sources have separate
licenses and runtime requirements. The upstream links below are acquisition
pointers, not an assertion that current upstream HEAD equals the executed
snapshot or that this release contains the frozen adapters.

| Method | Access pointer | Recorded call/settings and reference-count access |
|---|---|---|
| SNI-DPC | This repository | `SNIDPC().fit(X)`; one common fixed internal configuration |
| Gravity | GNaN [publication](https://doi.org/10.1016/j.patcog.2023.110063) and DSNGCAP [publication](https://doi.org/10.1016/j.asoc.2022.109647), including their linked implementation/supplementary records | Evaluated GNaN--DSNGCAP pipeline, not a separately bundled clustering package; its exact adapter is absent |
| Torque Clustering | [Official repository](https://github.com/JieYangBruce/TorqueClustering) | MATLAB Windows MEX, full Euclidean `DM`, `TORC(DM,0,0,0)`; first zero selects autonomous count mode |
| LDCC | [Publication record](https://doi.org/10.1016/j.patcog.2026.113065), with the authors' implementation/access information | Label-free adapter to published functions; `scale=10` for `n<500`, otherwise `100`; exact adapter not bundled |
| Gauging-delta | [Upstream repository](https://github.com/design-zeng/Gauging-delta) | Python `Perception(k=None)`; stop when no further merge occurs |
| DPC-MFP | [Upstream repository](https://github.com/Ancirled/DPC-MFP) | MATLAB/P-code `featureSelect.p`; common fixed defaults `k=42`, `l=19`; no reference count supplied |
| WANN-DPC | [Publication record](https://doi.org/10.1016/j.patcog.2025.111953) and its implementation/access records | Integer non-noise reference count supplied; Correction=2.5, NeighborLimit=30 with native data-dependent cap; no additional Min-Max; registered terminal-cycle safeguard |

Verify runtime compatibility, entry-point signatures and licenses after
acquisition. Keep MATLAB method paths isolated; one method's `normalize`
function must not shadow another's. The evaluated Torque MEX and DPC-MFP
P-code require their compatible MATLAB environment; an unrelated Python port
is not the same evaluated implementation. Exact comparator snapshot pins
and complete runnable configurations are not distributed in this repository.

## 4. Evaluation, repeats and resource records

NMI, Hungarian-matched ACC, ARI and Purity are computed after clustering.
For explicitly marked reference noise (`-1`), the primary clean-core scores
exclude those reference rows but all samples still enter clustering.
The reference count excludes the noise category; `c_out` is counted from
all predicted assignments. SNI-DPC assigns retained low-confidence points to
clusters and does not return an independent noise label. Full-sample noise
diagnostics, if used, have a distinct evaluation scope.

Each current method/dataset cell is one recorded call. Readable execution
paths do not use active algorithmic random sampling. The opaque Torque MEX
and DPC-MFP P-code had a separate ten-MATLAB-RNG-state check on one fixed
input; this was not a ten-repeat performance benchmark for every dataset.
The statistical scripts treat datasets, not RNG states, as pairs. WANN's
cycle safeguard is triggered on `solar_flare`; the other 38 predictions
remain pointwise unchanged under the same adapter. Its real-data resource
records come from independently timed single-worker calls, not a same-batch
seven-method rerun. The common configuration and source identities are in
[the evaluation record](seven-method-evaluation-config.json).

Recorded study environment: Windows 11 x64, Ryzen 9 7945HX (16 physical cores,
32 logical processors), 47.69 GB visible memory, Python 3.10.19 and MATLAB
R2026a Update 5. The Python environment is listed in
`requirements-experiments.txt`. This is study metadata, not a hardware
requirement or a claim that all archived records share one instantaneous load.

Runtime is the complete adapter call's wall-clock interval; peak memory is
the monitored process-tree RSS. The `page_blocks`/Gauging-delta record exceeded
the experimental script's configured **8-hour** limit and is retained as a
timeout (about 8.2 hours recorded elapsed time). Its quality/count fields are
empty, not zero. Original raw resource logs and per-task retry/configuration
history are not bundled. Budget settings from controlled sample-size,
dimensionality or noise tracks must not be inferred for these benchmark rows.

For what can be run directly from this repository, see
[Released materials](release-materials.md) and
[Current seven-method statistical protocol](seven-method-statistical-protocol.md).
The [six-method first-revision protocol](current-statistical-protocol.md)
and original scripts remain available with their own version identifiers.
