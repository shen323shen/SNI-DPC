# Data and comparator provenance

## Artificial benchmark

The formal benchmark contains 21 artificial datasets. Their source records and
the known provenance limitations are documented in the Supplementary Material.
This repository distributes the authoritative six-method result table, but does
not redistribute every inherited input file when the historical file identity,
version, or license cannot be established.

## Real-world study

The eight real-world cases are a bounded external-validity study. The exact
records used by the executed manifest are:

| Dataset role | Executed record or limitation |
|---|---|
| Image Segmentation | UCI data ID 50 |
| Soybean | OpenML data ID 42 |
| CPMP 2015 classification | OpenML data ID 41701 |
| Stock classification | OpenML data ID 841 |
| Tokyo1 | OpenML data ID 40705 |
| Solar Flare | OpenML data ID 40687 |
| Wine | inherited local copy; UCI data ID 109 is the upstream correspondence, not a byte-identity proof |
| Banknote | inherited local copy; UCI data ID 267 is the upstream correspondence, not a byte-identity proof |

OpenML data IDs identify concrete uploads. Same-name records, mirrors, and
current package versions must not be substituted for the executed record.
`xclara` is supplied by the R `cluster` package, while `hypercube` is generated
by `mlbench`; the historical package versions used for inherited files remain
unknown and are not silently replaced by current releases.

## Comparator code

The result table includes DPC-oracle-c, WANN-DPC, DSNGCAP, Gravity+DSNGCAP, and
KMeans++ (oracle-c). Only SNI-DPC source is distributed here. In particular,
DSNGCAP and Gravity+DSNGCAP are local evaluated combinations backed by group
publications; their source and license terms are not assumed to transfer to this
repository.

