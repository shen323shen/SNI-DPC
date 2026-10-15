# Data and comparator provenance

The current panel contains 29 artificial and ten real datasets.
[benchmark-input-manifest.json](benchmark-input-manifest.json) supplies
dataset-level identifiers, dimensions, counts and archived numeric CSV
SHA-256 values. [Acquisition and protocol](acquisition-and-protocol.md)
explains access and preparation. Inputs and prediction labels are not bundled.

## Artificial benchmark sources

| Dataset(s) | Source family or artifact |
|---|---|
| 2d-10c; triangle1 | Handl--Knowles generator family; Garcia-Piquer et al. (2015), [Information Sciences](https://doi.org/10.1016/j.ins.2015.04.041) |
| DS-577 | Su and Chou's fuzzy-clustering benchmark family (2005); historical original download no longer retained |
| R15 | Veenman et al., maximum variance cluster algorithm (2002) |
| banana | Immutable ClueMiner [banana1 ARFF artifact](https://github.com/deric/clueminer-fixtures/blob/2d694d56f66b53f6182ae8c2461e9df28b60c710/src/main/resources/data/bananas/banana.arff) |
| cluto-t4-8k; cluto-t5-8k; cluto-t8-8k | Karypis, CLUTO benchmark family (2002) |
| cure-t2-4k | Guha et al., CURE (1998), and the [type-2 generator](https://github.com/deric/handl-data-generators/tree/2a936ef96e6be180558de8609f370ef91a91c009) |
| diamond9 | Salvador and Chan, determining the number of clusters (2004) |
| dpb | Rodriguez and Laio, DPC Figure 2B (2014) |
| ds2c2sc13 | Faceli et al., benchmark partitions (2010) |
| flame | Fu and Medico, FLAME (2007) |
| fourty; longsquare; sizes4; smile1; smile3; twenty | Handl--Knowles / MOCK benchmark family |
| sizes2 | Handl and Knowles, TR-COMPSYSBIO-2004-02 and the MOCK family |
| hypercube; 2dnormals | R [mlbench](https://cran.r-project.org/package=mlbench) |
| xclara | R [cluster](https://cran.r-project.org/package=cluster) data object |
| zelnik5 | Zelnik-Manor and Perona, self-tuning spectral clustering (2004) |
| disk-6000n | [Immutable disk ARFF](https://github.com/deric/clustering-benchmark/blob/5a1e33ddf92e9ab6331f0c37c9589b4ff4818925/src/main/resources/datasets/artificial/disk-6000n.arff) and the Handl-data-generators source |
| bng_segment | [OpenML 130](https://www.openml.org/d/130), version 1; fixed 600-row processed sample |
| bng_mushroom | [OpenML 120](https://www.openml.org/d/120), version 1; fixed 600-row processed sample |
| aggregation | Gionis et al., [Clustering aggregation](https://doi.org/10.1145/1217299.1217303) (2007) |
| DS6 | [Immutable four-label MAT mirror](https://github.com/milaan9/Clustering-Datasets/blob/c09f0d17b9240900550ec46297d762d176e7908e/02.%20Synthetic/DS6.mat), 800 samples |

For the inherited artificial file collection, the manifest also links to an
immutable [Clustering Benchmark index](https://github.com/deric/clustering-benchmark/tree/5a1e33ddf92e9ab6331f0c37c9589b4ff4818925/src/main/resources/datasets/artificial).
A family/index correspondence is not an assertion that each newly downloaded
file is byte-identical to the archived numeric CSV.

Important file-level distinctions:

- `2d-10c`: the inherited ARFF declares ten classes; the evaluated CSV has
  nine nonempty reference labels and uses reference count nine.
- `banana`: the 4,811-by-2 two-class `banana1` artifact is traceable.
  Its original creator, generator and first publication remain unresolved.
  The mirror maintainer is not identified as the original dataset author.
- `DS6`: the evaluated version is 800-by-2 with four labels, not other
  same-name six-class or different-size variants. Original authorship and
  generator remain unresolved; the linked mirror identifies the artifact.
- `disk-6000n`: the source is consistent with
  `./disk -n 6000 -t 0 -b 0.2 -g 0.04`; its seed is not retained.
- `cure-t2-4k`: the type-2 generator defines 5% noise; its inherited seed
  is not retained.
- `2dnormals`: retained generation metadata is R `mlbench`,
  `set.seed(1234)`, `mlbench.2dnormals(1000)`.
- `hypercube`: 800 evaluated samples despite a retained 1,000-point header;
  historical package version is unverified. The historical `xclara`
  package version is also unverified.

## Real-data distributions

| Dataset | Persistent source | Version | File ID |
|---|---|---:|---:|
| banknote | [UCI 267](https://doi.org/10.24432/C55P57) | inherited local copy | not retained |
| cardiotocography | [OpenML 1466](https://www.openml.org/d/1466) | 1 | 1586231 |
| segmentation | [UCI 50](https://doi.org/10.24432/C5GP4N) | combined data/test files | not applicable |
| solar_flare | [OpenML 40687](https://www.openml.org/d/40687) | 4 | 4965264 |
| page_blocks | [OpenML 30](https://www.openml.org/d/30) | 1 | 30 |
| tokyo1 | [OpenML 40705](https://www.openml.org/d/40705) | 1 | 4965306 |
| soybean | [OpenML 42](https://www.openml.org/d/42) | 1 | 42 |
| cpmp_2015_classification | [OpenML 41701](https://www.openml.org/d/41701) | 1 | 21377443 |
| steel_plates_fault | [OpenML 40982](https://www.openml.org/d/40982) | 3 | 18151921 |
| stock | [OpenML 841](https://www.openml.org/d/841) | 2 | 53375 |

The banknote DOI is an upstream correspondence, not a byte-identity proof
for the inherited CSV. For segmentation the retained table appends the test
file after the training file. OpenML IDs/versions identify source uploads;
historical encoding, missing-value processing and row sampling can change
the numeric representation. Refer to the manifest for that representation's
shape and archived CSV hash, and to the acquisition guide for remaining
reconstruction requirements.

## Comparator code and historical version

The current six-method calls are SNI-DPC, Gravity, Torque Clustering, LDCC,
Gauging-delta and DPC-MFP. Only SNI-DPC source is distributed here.
[The comparator access table](acquisition-and-protocol.md#3-obtain-comparators-and-retain-the-documented-calls)
identifies upstream pointers, evaluated runtimes and known call settings.
Third-party redistribution rights are not inferred from public visibility.

DPC-oracle-c, WANN-DPC, DSNGCAP, Gravity+DSNGCAP and KMeans++ occur only in
the preserved historical v9.5 comparison and its script. They are not the
methods used by `reproduce_current_statistics.py`.
