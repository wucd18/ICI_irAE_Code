# Repository rules for code assistants

This repository stores workflow source code only, not research data or generated results.

- Do not modify model definitions, gene membership, filters, seeds, FDR families or biological sample selection during organization.
- Do not run historical all-in-one launchers, import scripts with top-level I/O, or regenerate approved paper figures to test a folder rearrangement.
- Keep data, fitted tables, final figures, workbooks, manuscripts and publication archives out of Git. Small script/config metadata maps are permitted.
- Respect `docs/script_index.csv`: `sha256` records original identity and `adapted_sha256` records delivered identity. Portability edits require a before/after diff and updated provenance, not an unrecorded replacement.
- Stage/execute in a separate workspace. Resolve original local imports and `__file__` roots explicitly; do not derive a workspace from arbitrary current working directories.
- Missing required external inputs must stop the stage, not cause synthetic or constant result generation.
- Never mark a stage as reproduced because syntax, source hashes or file presence passed. State each level of validation honestly.
- No remote GitHub or Figshare changes without explicit author approval. Preserve v1 and existing licences/history.
