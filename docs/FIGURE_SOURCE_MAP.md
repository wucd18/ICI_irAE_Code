# Figure and table source versions

| Version/topic | Real retained source | Confirmed from source | Still unverified |
|---|---|---|---|
| Recorded main figures | `plotting/recorded_closeout/plot_final.R`, `theme_final.R`, `csv_io.R` | Six recorded figures; 180 mm width; F6B contains a descriptive OLS line/band | Complete approved-version equivalence and external snapshot bindings |
| Later local main subset | `plotting/final_20260917/plot_final.R`, `theme_final.R`, `registered_original_csv_io.R` | Recovered from the final local build; width 170 mm; Figs 2/3/5/6; `mode=final` | Figs 1/4 final source closure, publication-version/layout identity; no replay |
| Recorded selected supplement | `plotting/recorded_closeout/plot_supplement_final.R` and `plot_organoid_supplement.R` | Distinct source-input/layout contracts | Final supplement mapping, not inferred from filename alone |
| Later S5/S7 | `plotting/supplement_20260917/prepare_plot_data.py`, `render_figures.py`, adapted `common.py` | Real CSV-to-display preparation and PDF-vector renderer | Inputs, fonts, exact final layouts, approved replay |
| Targeted extension/historical main | `plotting/targeted/*`, `plotting/historical_targeted/plot_v1.R` | Preserved distinct historical versions | Not automatically promoted to the final publication |
| Historical main/tables | `plotting/historical/make_manuscript_figures.R`, `plotting/tables/historical/build_final_synthesis_tables.py` | Numeric annotations/decisions now external inputs; original code retained locally | Exact annotations/archived report input binding |
| Recorded workbook | `plotting/tables/recorded/prepare_workbook.py`, `build_workbook.mjs` | True reader-table preparation and spreadsheet renderer | Public runtime portability and real cell/layout equivalence |
| Later workbook I/O | `plotting/tables/final_20260917/finalize_workbook.py`, adapted `common.py` | Real selective OOXML transfer with cell-value invariance checks | Exact external workbook/intermediate/edit-ledger binding and replay |

The two `plot_organoid_supplement.R` files are retained under different origins and stage families. They are never deduplicated by basename. All origin and adapted hashes appear in `script_index.csv`.

The required scientific architecture (study → receiver remodelling → frozen programme transfer → supportive TCR → organoid/TNFSF12 → final guardrail decision) is a scientific contract. Historical file/panel numbering is preserved as evidence; it was not rewritten to suggest that every archived publication version already matches that contract. Figure 6 remains the required final decision figure for any future coherent release. No manuscript, figure number or approved figure was edited in this task.
