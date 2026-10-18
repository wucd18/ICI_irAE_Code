# Scientific-source lineage and execution boundaries

## Patient results

The recovered recorded_upstream sources include patient-level and released-sample pseudobulk R analyses, fixed-programme derivation/validation, patient descriptive AUC Python, conditional R analysis, Oxford/Xenium Python and original organoid R/Python. Their original paths, hashes and available recorded script-index hashes are preserved in source_code_manifest.json.

Current publication plots read approved result tables. They do not rerun those scientific scripts. Recorded R scripts expect the original project configuration layout (including 00_scripts/project_config.R) and external inputs; WCD_PROJECT_ROOT, ICI_INPUT_ROOT and ICI_WORK_ROOT retain their original meanings. The supplied recorded configuration is provenance, not a claim that every upstream dependency has been bundled.

## All-ligand analysis

The approved source is CSBJ_ALL_LIGANDS_20261016T152734Z_r3. Actual DESeq2 code is increment_models.R; enrichment recovery is recover_enrichment.R and its Python orchestration/combination code.

Retained run history records that 257 gene-model/statistic files were reused unchanged from r2 after an enrichment-input UTF-8 BOM parsing failure. The r3 recovery repaired programme loading and used the preserved Wald rankings rather than fitting those models again. Historical worker completion and Windows runtime-exit diagnostics are retained in the original run archive. This history is not relabelled as a new end-to-end run.

The parameterized shared Python common.py needs:
ICI_ANALYSIS_RUN_ROOT (new writable destination), ICI_PROJECT_ROOT,
ICI_ORGANOID_SOURCE_ROOT, ICI_RECORDED_CODE_ROOT,
ICI_ANALYSIS_TASK_ROOT and ICI_FROZEN_MEMBERSHIP.
R entrypoints accept the analysis output root and shard arguments and read analysis_config.json.
Standalone recovery/diagnostic/orchestration scripts retain some original paths and require deliberate review before any future execution.

## Latest directional/background sensitivity

The source is CSBJ_DIRECTIONAL_BACKGROUND_20261016T183734Z.
prepare_context.py supplies matched-control baselines; build_libraries.py supplies matched gene libraries; main_scores.py preserves original, background-relative and non-rank responses; pool_sensitivity.py and composition_sensitivity.py supply the established pool/common-support summaries. All are actual retained sources.

Its shared common.py needs:
ICI_ANALYSIS_RUN_ROOT (new writable destination), ICI_ALL_LIGAND_RUN
(approved saved analysis), ICI_ORGANOID_SOURCE_ROOT and
ICI_SENSITIVITY_TASK_FILE.
The analysis_amendment.md records the existing sensitivity specification.
Original standalone document/diagnostic wrappers may still contain historical machine paths.

## What was executed now

Publication plotting/assembly, table-result formatting, Word/PDF production and source syntax checks only. No edgeR, DESeq2 or fgsea model/permutation run was started, no scores were newly defined, and no gene, donor, ligand, eligibility or correction family was changed.

R syntax was checked using Rscript --vanilla parse(file=...), without sourcing or evaluating scientific code. Existing Windows locale-startup warnings were recorded locally; all nine R sources parsed successfully. Python and Node syntax checks are code checks, not biological reproduction.

External raw matrices, intermediate results, random libraries and model objects remain in the original project. This package supplies the actual current production source and its scientific lineage; independent raw-data end-to-end portability remains untested.

## Source-label uncertainty

IFNG and TNFA remain the exact source keys and map to IFN-γ and TNF in prose. The existing IL22/IL11 source-label/dose mapping uncertainty is preserved rather than guessed; the relevant table fields remain NA. Authors should resolve that original metadata identity before claiming a specific mapped dose.

