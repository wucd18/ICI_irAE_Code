#!/usr/bin/env python3
"""Create a checksum-backed manifest for every file in 01_data."""

from __future__ import annotations

import argparse
import hashlib
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


SUPPLEMENT_URLS = {
    "GSE228597_combined_tissue_data.h5ad.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE228nnn/GSE228597/suppl/GSE228597_combined_tissue_data.h5ad.gz",
    "GSE253720_Biopsy_RNA_Final.h5ad.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253720/suppl/GSE253720_Biopsy_RNA_Final.h5ad.gz",
    "GSE210037_TG_RR_counts_raw.tsv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210037/suppl/GSE210037_TG_RR_counts_raw.tsv.gz",
    "GSE210037_TG_RR_sample_metadata.tsv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE210nnn/GSE210037/suppl/GSE210037_TG_RR_sample_metadata.tsv.gz",
    "GSE277136_BLFplusNM.h5ad": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277136/suppl/GSE277136_BLFplusNM.h5ad",
    "GSE277136_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277136/suppl/GSE277136_RAW.tar",
    "GSE273720_sample_list.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE273nnn/GSE273720/suppl/GSE273720_sample_list.csv.gz",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(16 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def infer_url(name: str) -> tuple[str, str]:
    if name in SUPPLEMENT_URLS:
        return SUPPLEMENT_URLS[name], "downloaded"
    if name.endswith("_family.soft.gz"):
        accession = name.split("_", 1)[0]
        series_bucket = accession[:6] + "nnn"
        return (
            f"https://ftp.ncbi.nlm.nih.gov/geo/series/{series_bucket}/{accession}/soft/{name}",
            "downloaded",
        )
    if name.endswith(".h5ad") and f"{name}.gz" in SUPPLEMENT_URLS:
        return f"derived_by_gzip_decompression_from:{name}.gz", "derived"
    return "local_or_unmapped", "local_or_unmapped"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    args = parser.parse_args()
    project = args.project.resolve()
    data_dir = project / "01_data"
    table_path = project / "tables_for_article" / "Table_S10_download_and_file_manifest.tsv"

    rows = []
    files = sorted(x for x in data_dir.iterdir() if x.is_file())
    for index, path in enumerate(files, start=1):
        url, provenance = infer_url(path.name)
        print(f"{index}/{len(files)} hashing {path.name}", flush=True)
        stat = path.stat()
        rows.append(
            {
                "file_name": path.name,
                "bytes": stat.st_size,
                "sha256": sha256(path),
                "provenance": provenance,
                "source": url,
                "file_modified_utc": datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
                "manifest_checked_utc": datetime.now(tz=timezone.utc).isoformat(),
            }
        )
    pd.DataFrame(rows).to_csv(table_path, sep="\t", index=False)
    print(f"Wrote {table_path} with {len(rows)} files")


if __name__ == "__main__":
    main()
