# Implementation equivalence

The public estimator imports `src/sni_dpc/_core_v10.py`, generated from the
preserved v9.5 core and the six revised functions in
`frozen/SNI_DPC_v10_0_65.py`. Its composite-density term assigns a fixed
weight of 0.65 to the structural-field view and 0.35/3 to each of the other
three views. The SNI-DPC screened structural field force is distinct from the
Gravity comparison method's field.

Run `python scripts/build_v10_core.py --check` to verify generator inputs and
the checked-in core. The reference temporarily patches functions in the
historical core during its call; use the public estimator for concurrent
applications. The generated module binds functions locally without patching
the historical module. Tests compare pointwise labels, predicted counts,
natural-neighbor order and low-confidence indices on deterministic inputs;
this does not claim a full 39-dataset benchmark rerun.

## Algorithm-to-code map

| Manuscript stage | Frozen implementation functions |
|---|---|
| Extended natural-neighbor search | `natural_neighbor_search`, `_query_knn_excluding_self` |
| Screened structural field | `build_mng`, `compute_reliability_from_nb`, `compute_local_scale_from_gamma`, `compute_structure_similarity`, `compute_gravity`, `compute_composite_density` |
| Low-confidence filtering | `detect_low_density_points_by_dsngcap_density`, `detect_low_density_points_by_nb_adaptive` |
| Clean cluster assembly | `cluster_clean_data`, `build_component_mng_by_strong_edges`, `select_centers_by_raw_cut_budget`, `layer1_core_expansion`, `layer2_prototype_expansion`, `layer3_gravity_boundary_assignment` |
| Output and reassignment | `mng_dpc`, `redistribute_low_density_points` |

The two low-confidence branches are combined by set union before clean-graph
reconstruction. This operation is not center merging. Center-budget estimation,
center scoring, and center detection occur later in the clean clustering stage.

The unchanged historical `frozen/SNI_DPC_v9_5_publication.py` and 21-dataset
table remain for traceability. Hashes for both generations are in
`checksums.sha256`; historical oracle-count baselines do not define the
current comparison.
