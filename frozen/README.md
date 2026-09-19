# Checked reference implementations

`SNI_DPC_v10_0_65.py` records the fixed-weight v10_0.65 reference used for
public-API equivalence checks. It imports the historical core and temporarily
overrides selected functions during each call; use the public `SNIDPC`
estimator for applications, especially when calls can overlap. Its hash is
recorded in `../checksums.sha256`.

`SNI_DPC_v9_5_publication.py` is the unchanged historical submission
snapshot. Its hash is also recorded. Do not reformat either checked reference;
see `../docs/implementation-equivalence.md` for their relationship.
