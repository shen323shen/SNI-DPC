# Implementation equivalence

The file `frozen/SNI_DPC_v9_5_publication.py` is the publication copy of the
implementation used to produce the submitted SNI-DPC results. Its executable
body was previously verified line-for-line against the frozen experimental
source. The extra publication header does not change computation.

The immutable SHA-256 values are recorded in `checksums.sha256`:

- frozen publication implementation:
  `BBD3FD61314325654882C4EDF158B70A9F4CE11E19E7FD64298F30A99EFE4671`;
- formal 21-dataset, six-method result table:
  `2D6518035E79D9DA5397D05B64D9DDF0C3E392F0AF8BE812F9AF6097820F44EF`.

## Manuscript-to-code map

| Manuscript stage | Frozen implementation functions |
|---|---|
| Natural-neighbor search | `natural_neighbor_search`, `_query_knn_excluding_self` |
| Screened interaction field | `build_mng`, `compute_reliability_from_nb`, `compute_local_scale_from_gamma`, `compute_structure_similarity`, `compute_gravity`, `compute_composite_density` |
| Low-confidence filtering | `detect_low_density_points_by_dsngcap_density`, `detect_low_density_points_by_nb_adaptive` |
| Clean cluster assembly | `cluster_clean_data`, `build_component_mng_by_strong_edges`, `select_centers_by_raw_cut_budget`, `layer1_core_expansion`, `layer2_prototype_expansion`, `layer3_gravity_boundary_assignment` |
| End-to-end clustering | `mng_dpc`, `redistribute_low_density_points` |

The two low-confidence branches are combined by set union before clean-graph
reconstruction. This operation is not center merging. Center-budget estimation,
center scoring, and center detection occur later in the clean clustering stage.

The cleaned package under `src/sni_dpc/` is tested against this frozen reference.
The frozen file is intentionally excluded from automatic formatting.
