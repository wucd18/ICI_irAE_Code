# Current figure and table source map

All production paths below are within production/final_20261018/code. The display input root is external and is recorded in local_input_locations.json.

| Publication output | Actual current generator | Inputs and interpretation |
|---|---|---|
| Figure 1 | plot_main.py, assemble_main.py | Fixed programme definitions, patient coverage and patient-level scores from source_data/old |
| Figure 2 | plot_main.py, assemble_main.py | Patient enrichment, descriptive AUC and Oxford results from source_data/old |
| Figure 3 | redraw_figure3.py, assemble_main.py | Complete planned comparison coverage, three representations and donor effects; approved v1_panel_data used where specified by the recovered code |
| Figures 4–6 | redraw_main_4_to_6.py, assemble_main.py | source_data/new donor changes, member effects, composition and reference overlap/correlation tables |
| Figure S1 | plot_supplement_final.py, assemble_supplement_final.py | External and discovery programmes, AUC and individual patient scores |
| Figure S2 | same | Existing IFN-adjusted associations and observed support; residual c–g above depleted h–l |
| Figure S3 | same | Oxford four-programme assessment and Xenium descriptive epithelial scores |
| Figure S4 | same | One complete non-rank continuous ligand-response atlas across the four families |
| Figure S5 | same | Direction, matched-background and experimental-pool sensitivity summaries |
| Figure S6 | same | All original/additional Both combinations and score-dependent donor changes |
| Figure S7 | same | Recovered composition, fixed-state-weight responses and remaining-family reference scatterplots |
| Tables S1–S9 | prepare_tables.py, build_workbook.mjs | Complete defined publication sections from saved old/new results plus pinned IBD memberships |
| Manuscript and Table 1 | revise_manuscript.py | Exact reviewed ZIP baseline plus bounded text/reference edits and refreshed figure previews |
| Manuscript_with_Figures.pdf | export_documents.ps1, assemble_manuscript_pdf.py | Actual Word export plus six native vector figure pages |
| Cover_Letter.docx | preserved baseline file | Authors, reviewers and existing 16 October date retained |

The current workbook is byte-identical to the reviewed workbook. Its genuine previous generator is included and its newly parameterized table-source preparation was executed with identical content. No spreadsheet data were recomputed.

theme_01.py and style_01.json control stable family colours, donor shapes and local fonts; figure_layout.json controls panel layout. render_main_rasters.py supplies Word previews. The two contrast fixes choose black or white label text on unchanged S5e/S6a classification colours. All NE cells, donor points, original ranges and result values remain.

The earlier 25-page supplement/nine-workbook layout is not the current production route. Older plot and prose scripts under the analysis folders are retained as scientific-run source history; use only production/final_20261018 for this publication layout.

