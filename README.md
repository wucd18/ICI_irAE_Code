# ICI-irAE analysis code

**Code-only upload candidate — PORTABILITY_PENDING. This is not a validated end-to-end reproduction package.**

## Study
These sources examine shared inflammatory pressure, context-dependent receiver-cell remodelling, disease-depleted epithelial programmes, and apparent transcriptional-restoration signals against tumour-response and off-target guardrails. Patients or released samples are the biological units. TCR is supportive evidence, not proof of an epithelial mechanism. AUC assesses repeatability, not diagnostic performance; absence of detected tumour harm does not establish pharmacologic safety.

## Code only
This repository contains analysis, data preparation, plotting and table-generation programs, configuration metadata and instructions. It contains no actual research results, finished figures, result workbooks, expression matrices, manuscript or data ZIP. `plotting/` contains programs, not pictures; `plotting/tables/` contains programs, not result tables. Results arise only from separately authorized runs with real inputs.

## Directories
| Directory | Contents |
|---|---|
| data_prep | Acquisition, patients, organoid exports/counts, conditional, spatial, and extracted display-input preparation |
| analysis | Patient-level models, programmes, organoid, conditional, TCR, clinical and spatial; historical selection is separately labelled |
| plotting | Recorded and historical renderers, later recovered figure versions, and table I/O |
| config | Explicit path example and deliberately incomplete external-input manifest example |
| environments | Source-derived dependency lists and the original presentation R session record |
| workflow | Read-only planning, input hash/schema checks, and fresh external code staging |
| tests | Engineering fixtures/static checks; archived scientific verifiers are separately labelled |
| docs | Source identities, dependencies, external-input requirements, figure map and limitations |

## Data versus code
`INPUT_ROOT` and `ARCHIVE_ROOT` are read-only external locations. Frozen gene memberships, patient mappings, saved estimates, historical decisions and display annotations are external inputs, even when small. The author-supplied [Figshare locator](https://doi.org/10.6084/m9.figshare.33872563) is a locator only: the exact remote version and member identities have not been verified in this build. Local fingerprints do not certify remote completeness.

## Installation
No packages were installed or upgraded. `environments/python_packages_from_imports.txt` and `r_packages_from_source.txt` are unpinned dependency candidates derived from code, not a freeze or lock. After approval, create an isolated Python environment and install the reviewed requirements there. R environments require separate CRAN/Bioconductor setup compatible with each recorded version. See [environment limitations](environments/README.md); do not silently unify the original R 4.2.2 launcher reference and the targeted/presentation R 4.4.1 records.

## Configuration
Copy `config/paths.example.json` to an ignored local config and supply absolute `CODE_ROOT`, `INPUT_ROOT`, `ARCHIVE_ROOT`, `WORK_ROOT` paths. Code and input roots must not overlap the work root. Executables are explicit; the tool does not discover private environments. Optional `OUTPUT_ROOT` must be beneath `WORK_ROOT`. Put a real versioned input manifest outside Git; the empty example intentionally fails verification.

## Concrete commands
From this repository, using your configured Python executable:
```text
python -B tests/check_repository.py --r-parse --rscript /absolute/path/to/Rscript
python -B tests/test_portability.py
python -B workflow/stage_workspace.py --help
python -B workflow/stage_workspace.py --config /absolute/path/paths.local.json --plan
python -B workflow/stage_workspace.py --config /absolute/path/paths.local.json --verify-inputs /absolute/path/verified_inputs.json --plan
python -B workflow/stage_workspace.py --config /absolute/path/paths.local.json --write
```
`--write` stages code only in `WORK_ROOT/runs/<unique_id>/`; it does not copy data, run models, render figures, or install software. Plans and help are read-only. New staging refuses existing run directories. Do not execute top-level legacy programs in the classified source tree. [RUN_ORDER.md](docs/RUN_ORDER.md) provides stage-specific command contracts and the remaining gates. There is no advertised one-command raw-to-final runner.

## Figure and script mapping
[FIGURE_SOURCE_MAP.md](docs/FIGURE_SOURCE_MAP.md) distinguishes the recorded 180 mm six-figure renderer from later 170 mm subset sources, supplementary S5/S7 sources and workbook versions. Recovered code identity does not prove equivalence to every approved figure or table. No approved figure was regenerated or replaced.

## Reproduction levels
Static syntax, source identity and engineering fixture checks were executed locally. Input fingerprints were checked for selected local candidates; comprehensive stage schemas and remote archive mapping remain incomplete. Model execution, numerical reproduction, full raw-to-final replay, and figure/workbook equivalence were **NOT_RUN**. A static PASS must never be relabelled as scientific reproduction.

## Known limits
See [PORTABILITY.md](docs/PORTABILITY.md) for exact gaps: mixed legacy read/write layouts, CSV/TSV bridges, versioned input/member binding, externalized historical annotations/reference membership, missing local dependencies, the proprietary workbook runtime dependency, and incomplete final figure/layout correspondence. Fixed thresholds, seeds, model definitions, candidate gates and biological units were not tuned. Historical selection and donor re-evaluation are not new independent experiments. Historical machine identifiers such as `FUNCTION_RESCUE_CORE` remain identifiers; use “epithelial reciprocity core” in interpretation. SORL1 remains exploratory.

## Licence and contact
No licence or software-contributor list was invented. [LICENSING.md](docs/LICENSING.md) records the ownership/permission gap. The owner must retain any legitimate existing licence and confirm redistribution rights before publication. No CITATION.cff was fabricated. The supplied repository locator is [wucd18/ICI_irAE_Code](https://github.com/wucd18/ICI_irAE_Code); no remote repository, tag, release or manuscript was changed. No contact identity was inferred.
