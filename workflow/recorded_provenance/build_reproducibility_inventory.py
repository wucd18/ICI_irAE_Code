#!/usr/bin/env python3
"""Create checksum inventories for analysis scripts and final submission artifacts."""

from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rows_for(paths, root):
    rows = []
    for path in sorted(paths):
        stat = path.stat()
        rows.append({
            "relative_path": path.relative_to(root).as_posix(),
            "bytes": stat.st_size,
            "sha256": sha256(path),
            "modified_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
        })
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("project_dir", type=Path)
    args = parser.parse_args()
    root = args.project_dir.resolve()
    tables = root / "tables_for_article"

    scripts = [p for p in (root / "00_scripts").iterdir() if p.is_file()]
    script_df = pd.DataFrame(rows_for(scripts, root))
    script_df.insert(1, "language", script_df["relative_path"].str.rsplit(".", n=1).str[-1])
    script_df.to_csv(tables / "Table_S62_reproducibility_script_inventory.tsv", sep="\t", index=False)

    artifacts = []
    for dirname in ["figures_for_article", "tables_for_article", "manuscript"]:
        artifacts.extend(
            p
            for p in (root / dirname).iterdir()
            if p.is_file() and p.name != "Table_S63_final_artifact_manifest.tsv"
        )
    artifact_df = pd.DataFrame(rows_for(artifacts, root))
    artifact_df.to_csv(tables / "Table_S63_final_artifact_manifest.tsv", sep="\t", index=False)
    print(f"Inventoried {len(script_df)} scripts and {len(artifact_df)} final artifacts")


if __name__ == "__main__":
    main()
