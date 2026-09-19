# Data and comparator provenance

## Artificial benchmark

The current six-method panel has 29 artificial datasets: the earlier 21 plus
eight extensions (`sizes2`, `disk-6000n`, `cluto-t8-8k`, `bng_segment`,
`bng_mushroom`, `2dnormals`, `aggregation`, `DS6`). The published CSV contains
one row per method and dataset. Source records, file identity and any
historical provenance limitations must be checked against the executed
dataset manifest; this repository does not claim that a same-name download
is byte-identical to the evaluated local file. Inputs without verified
redistribution terms are not bundled.

## Real-world study

The current comparison includes ten real-world cases. Dataset identifiers
recorded in the accompanying study are shown below; a matching public ID
alone does not prove byte identity with a locally executed copy.

| Dataset role | Executed record or limitation |
|---|---|
| banknote | inherited local copy; UCI data ID 267 is an upstream correspondence, not a byte-identity proof |
| cardiotocography | OpenML data ID 1466 |
| segmentation | UCI Image Segmentation, data ID 50 |
| solar_flare | OpenML data ID 40687 |
| page_blocks | OpenML data ID 30 |
| tokyo1 | OpenML data ID 40705 |
| soybean | OpenML data ID 42 |
| cpmp_2015_classification | OpenML data ID 41701 |
| steel_plates_fault | OpenML data ID 40982 |
| stock | OpenML data ID 841 |

OpenML IDs identify concrete uploads. Same-name records, mirrors and current
package versions cannot automatically replace executed input files. For the
historical artificial inputs, `xclara` is supplied by R `cluster` and
`hypercube` by `mlbench`; inherited package versions remain unverified.

## Comparator code

The current tables compare SNI-DPC, Gravity, Torque Clustering, LDCC,
Gauging-delta and DPC-MFP. Only the SNI-DPC source is distributed here;
third-party implementation licensing is not inferred from the existence of
public repositories. DPC-oracle-c, WANN-DPC, DSNGCAP, Gravity+DSNGCAP and
KMeans++ occur solely in the preserved historical 21-dataset table; they are
not current comparison methods.
