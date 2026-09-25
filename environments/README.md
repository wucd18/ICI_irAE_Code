# Environment evidence

Requirements were extracted from real Python imports and R package references. They are unpinned. `dependency_references.csv` retains the calling file and observed local availability. Its installed versions describe the configured local interpreter only; no historical Python freeze or R lockfile is claimed.

The original launcher in code.zip references R 4.2.2. The supplied `recorded/presentation_R_sessionInfo.txt` records R 4.4.1 for presentation. Targeted source records R 4.4.1-era dependencies. These have not been unified or recreated. Local syntax parsing used the executable in the parent project's `config/local_paths.ps1`, with `--vanilla`; the detailed actual versions and library paths are in LOCAL_REVIEW_ONLY, outside this public tree.

At inspection, the configured Python environment lacked beautifulsoup4, lxml, openpyxl and reportlab; the queried R library paths lacked metafor. R package presence alone is not loadability or numerical equivalence. The Windows environment probe is retained historical source and may perform network requests; it was not executed. New staging never reads/restores a private environment JSON.

No installation, upgrade, environment replacement or network access was performed. Example setup commands for an independently approved setup step are `python -m venv <external-env>` followed by that environment's `python -m pip install -r environments/python_packages_from_imports.txt`. These unpinned commands are not a reproduction guarantee. R/Bioconductor setup must be planned for the applicable historical version; do not install/upgrade during an analysis run.

The workbook renderer still uses the artifact-tool API. `ICI_ARTIFACT_TOOL_MODULE` must be an explicit legitimately installed module file; its public installation availability was not asserted. A different public renderer requires validated cell-value/order/missingness and layout comparisons before claiming equivalence. Arial font files for the recovered supplementary renderer must be provided locally under appropriate rights using explicit font paths; fonts are not redistributed.
