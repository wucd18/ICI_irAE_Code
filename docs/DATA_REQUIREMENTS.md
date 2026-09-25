# External input contracts

| Stage | Required real input families |
|---|---|
| Patient inference | Pseudobulk counts/genes/metadata and actual donor/condition/cell-family mappings; reviewed raw input metadata |
| Frozen programme transfer | `Table_S21_frozen_module_definitions.tsv`, `Table_S13_GSE206300_cluster_family_mapping.tsv`, original score inputs; membership stays outside Git |
| Historical selected sets | `Table_S25_ICB_guardrail_feature_definitions.tsv` or a provenance-verified corresponding Source table; use `up_genes`/`down_genes` for the historically selected signatures/core, never infer members from a figure |
| Organoid | `organoid_donor_extension/source_data/pseudobulk_counts.mtx`, `pseudobulk_metadata.tsv`, `pseudobulk_genes.tsv`; matching locked plan/QC; original fixed definitions |
| Conditional | Real patient scores/counts, frozen definitions/family mapping, locked conditional plan; the saved plan is not a success switch |
| Spatial | Released count layers, patient/compartment/section metadata, coordinates where released, frozen genes and source manifests |
| Display | Distinct verified snapshots, current/final panel tables, selected-input manifest, exact recorded organoid result directory |
| Historical annotation/reference | `ICI_HISTORICAL_DISPLAY_LABELS`, `ICI_HISTORICAL_ANNOTATIONS`, `ICI_RECEPTOR_MAP`; schema keys are in `external_annotation_requirements.csv` |
| Workbook | Recorded reader tables and table definitions; later workbook requires original XLSX plus authorized cell-change/intermediate records |

No official data manifest is invented here. `config/input_manifest.example.json` is intentionally empty and fails. Supply a version plus nonempty `files` entries containing `id`, `root` (INPUT_ROOT/ARCHIVE_ROOT), relative `path`, `size`, `sha256`, accession/version/source, role and stage. Tabular entries can declare `columns`, `unique_keys` and delimiter. The verifier rejects missing, changed, duplicated or escaping inputs. Comprehensive per-stage schema/patient cross-file validation is still required.

The local review contains byte fingerprints of found candidates, not scientific data copies. A matching schema does not prove version equivalence. Root DOI accessibility and version/member hashes were not checked online. Do not silently use another cache, prior-run result, empty table, synthetic table or a digitized figure to fill a missing formal input.
