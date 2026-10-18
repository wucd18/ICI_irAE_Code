# Running the current publication generators

## External bindings

Copy `local_bindings.json` to a private working configuration and adjust real paths. The file records the local paths used for this delivery; it contains no credentials. Inputs are external, read-only project assets:

| Binding | Required content |
|---|---|
| ICI_DISPLAY_INPUT_ROOT | Approved display_revision_v2, including source_data/old, source_data/new and archive/v1_panel_data |
| ICI_SUBMISSION_BASELINE | Exact submission_final extracted from ICI_SUBMISSION_FINAL_20261018.zip |
| ICI_IBD_MEMBERSHIP | Pinned IBD_membership_definitions.tsv from the approved all-ligand run |
| plot_python | Existing scientific Python with NumPy, pandas and Matplotlib |
| artifact_python | Existing document Python with pandas, pypdf, lxml and other imports used by the scripts |
| node / artifact_tool_module | Existing Node and the actual bundled @oai/artifact-tool module file |
| powershell | PowerShell able to invoke installed Microsoft Word COM |
| pdftoppm | Existing Poppler executable |
| ICI_ARIAL_REGULAR / ICI_ARIAL_BOLD (optional) | Local Arial TTF files; default Windows font paths are retained |

The baseline ZIP SHA-256 is
`d0f1899a2b6bd805d799baa81051d03a82171c81e5f1d26814d3f18bf35233b5`.
Use this baseline, not the already revised manuscript or a historical JTM document. The narrative editor is a deterministic closeout patch with paragraph identities tied to that baseline.

## Presentation command

From this code directory, use the Python executable specified by your bindings:

```powershell
& '<artifact_python>' production/final_20261018/rebuild_publication.py --bindings '<private-bindings.json>' --out '<new-output-directory>' --stage all
```

The output directory must not exist. It must be outside all input trees. No analysis model is invoked.

The wrapper creates a new submission_final and separate local working directories. Only submission_final belongs in the journal ZIP. Keep inputs, panels, panel_data, previews and qa outside that ZIP.

The available stages are:
- figures: rebuild current figures and supplement from saved numeric tables;
- tables: format all nine complete result tables;
- documents: apply the closeout text patch and use the baseline figures for a document-only build;
- all: rebuild figures and tables, then insert those current figures and export the documents.

Use all for a complete new publication build. The document-only mode intentionally does not redraw figures. Actual corresponding production stages were tested separately during closeout; no single-command fresh-machine scientific run is asserted.

Microsoft Word is used for line layout and pagination. assemble_manuscript_pdf.py retains the exported text pages and replaces the six figure pages with native vector PDFs. The Word copy contains the raster previews necessary for Word display, while the delivered figure PDFs retain vector graphics and embedded text.

The spreadsheet generator is the actual executed build_workbook.mjs. Its dependency is the preinstalled Codex bundled artifact-tool (recorded version in environment.json), not a declared public npm package. No install command is supplied. On a machine without that module or Microsoft Word, the relevant original production route remains an explicit environment dependency.

## Scientific sources

Do not run the analysis entrypoints for presentation changes. They can fit models or rewrite analysis output. A future authorised scientific run needs its own empty output root, verified external data, and the recorded R environment. See ANALYSIS_SCOPE.md for the source lineage and extra bindings. No historical run directory should be used as a writable destination.

