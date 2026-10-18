# Mapping recurrent epithelial programmes of checkpoint inhibitor colitis to human organoid ligand responses

Local final-code delivery, 18 October 2026.

The current publication contains six main figures, seven supplementary figures on 19 pages, and nine editable result worksheets in one Supplementary_Tables.xlsx. This directory accompanies Submission_Text_Final.zip; it is not an online release.

## What to use

- `production/final_20261018/`: the actual current figure, supplement, table, Word and PDF generators, with explicit local input/output bindings.
- `analysis/all_ligand/`: the retained all-ligand Python and R sources associated with the approved 16 October results.
- `analysis/directional_background/`: the latest retained matched-background, non-rank, pool, composition and reference-comparison Python sources.
- `recorded_upstream/`: relevant patient, programme, conditional, spatial and original organoid sources recovered from the existing recorded code tree.
- `docs/source_code_manifest.json`: original locations/hashes and the packaged adaptation of each recovered file.

The new `rebuild_publication.py` wrapper binds paths and invokes the recovered generators. It is explicitly new glue code, not a recovered historical analysis script. Narrative edits in `revise_manuscript.py` were created and executed for this closeout against the exact reviewed manuscript.

Read [RUNNING.md](docs/RUNNING.md), [FIGURE_SOURCE_MAP.md](docs/FIGURE_SOURCE_MAP.md), and [ANALYSIS_SCOPE.md](docs/ANALYSIS_SCOPE.md) before use.

## Scope and executed checks

The final production sources were used to rebuild all six figures and the 19-page supplement in a separate local output directory. Their plot-source values, PDF text and page geometry matched the delivery. The nine-table source generator produced identical table content. Actual Microsoft Word export produced the revised manuscript PDF.

Scientific model scripts were recovered and parsed, not rerun. Current visual and numerical correspondence checks are not raw-data end-to-end reproduction. Detailed QA, logs, rendered pages, prior files and data remain in the local project archive and are intentionally absent from this code ZIP.

No source file required for current publication production was missing. External saved results, frozen membership definitions and the exact reviewed document baseline are required; their local bindings are documented. The bundled spreadsheet module and Word runtime are also required for the current original output route.

## Publication and licensing

The existing project repository is https://github.com/wucd18/ICI_irAE_Code. This local delivery has not been pushed, tagged or published, and creates no release or DOI. It does not replace the public root README automatically. No new licence is granted by this packaging step; retain applicable source rights and obtain the authors' licence decision before a public code release. Fonts, raw data, model objects and publication outputs are not redistributed here.

