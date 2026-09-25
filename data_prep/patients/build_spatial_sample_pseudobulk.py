#!/usr/bin/env python3
"""Stream a spot-by-gene text matrix into patient/sample-level spatial pseudobulk.

The GEO file is gene-by-spot and too wide for a dense data-frame import. This
script reads one gene at a time, verifies every barcode against the supplied
metadata, and aggregates raw counts by biological sample without retaining a
dense spot-level expression matrix.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    args = parser.parse_args()
    project = args.project.resolve()

    data_dir = project / "01_data"
    results_dir = project / "02_results"
    counts_path = data_dir / "GSE210037_TG_RR_counts_raw.tsv.gz"
    metadata_path = data_dir / "GSE210037_TG_RR_sample_metadata.tsv.gz"

    metadata = pd.read_csv(metadata_path, sep="\t")
    required = {"barcode", "sample", "tissue", "nCount_Spatial", "nFeature_Spatial"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Missing spatial metadata columns: {sorted(missing)}")
    if metadata["barcode"].duplicated().any():
        raise ValueError("Spatial metadata contains duplicate barcodes")

    with gzip.open(counts_path, "rt") as handle:
        matrix_barcodes = handle.readline().rstrip("\r\n").split("\t")
        matrix_barcodes = [x.replace(".", "-", 1) for x in matrix_barcodes]
        metadata_by_barcode = metadata.set_index("barcode", drop=False)
        missing_barcodes = sorted(set(matrix_barcodes).difference(metadata_by_barcode.index))
        extra_barcodes = sorted(set(metadata_by_barcode.index).difference(matrix_barcodes))
        if missing_barcodes or extra_barcodes:
            raise ValueError(
                f"Barcode mismatch: matrix-only={len(missing_barcodes)}, "
                f"metadata-only={len(extra_barcodes)}"
            )

        ordered = metadata_by_barcode.loc[matrix_barcodes].reset_index(drop=True)
        samples = sorted(ordered["sample"].astype(str).unique())
        sample_index = {sample: i for i, sample in enumerate(samples)}
        group_codes = ordered["sample"].astype(str).map(sample_index).to_numpy(dtype=np.int64)

        genes: list[str] = []
        gene_sums: list[np.ndarray] = []
        for line_number, line in enumerate(handle, start=2):
            gene, separator, payload = line.rstrip("\r\n").partition("\t")
            if not separator:
                raise ValueError(f"Malformed counts row at line {line_number}")
            values = np.fromstring(payload, dtype=np.int64, sep="\t")
            if values.size != len(matrix_barcodes):
                raise ValueError(
                    f"Count length mismatch for {gene}: {values.size} != {len(matrix_barcodes)}"
                )
            genes.append(gene)
            gene_sums.append(
                np.bincount(group_codes, weights=values, minlength=len(samples)).astype(np.int64)
            )
            if len(genes) % 5000 == 0:
                print(f"aggregated {len(genes):,} genes", flush=True)

    gene_by_sample = np.vstack(gene_sums)
    sample_by_gene = sparse.csr_matrix(gene_by_sample.T)
    matrix_path = results_dir / "GSE210037_spatial_pseudobulk_counts.mtx.gz"
    with gzip.open(matrix_path, "wb") as output:
        mmwrite(output, sample_by_gene, field="integer", symmetry="general")

    gene_table = pd.DataFrame({"gene": genes, "gene_index_1based": np.arange(1, len(genes) + 1)})
    gene_table.to_csv(results_dir / "GSE210037_spatial_pseudobulk_genes.tsv", sep="\t", index=False)

    sample_metadata = (
        ordered.groupby(["sample", "tissue"], observed=True)
        .agg(
            n_spots=("barcode", "size"),
            metadata_total_counts=("nCount_Spatial", "sum"),
            median_features_per_spot=("nFeature_Spatial", "median"),
        )
        .reset_index()
        .set_index("sample")
        .loc[samples]
        .reset_index()
    )
    sample_metadata["matrix_total_counts"] = np.asarray(sample_by_gene.sum(axis=1)).ravel()
    sample_metadata["count_total_matches_metadata"] = (
        sample_metadata["matrix_total_counts"].astype(np.int64)
        == sample_metadata["metadata_total_counts"].astype(np.int64)
    )
    sample_metadata.insert(0, "row_index_1based", np.arange(1, len(sample_metadata) + 1))
    sample_metadata.to_csv(
        results_dir / "GSE210037_spatial_pseudobulk_metadata.tsv", sep="\t", index=False
    )

    audit = {
        "dataset": "GSE210037",
        "matrix_orientation": "samples_by_genes",
        "n_samples": int(sample_by_gene.shape[0]),
        "n_genes": int(sample_by_gene.shape[1]),
        "n_spots": int(len(matrix_barcodes)),
        "nnz": int(sample_by_gene.nnz),
        "barcodes_exactly_matched_after_first_dot_to_hyphen": True,
        "all_sample_totals_match_spot_metadata": bool(
            sample_metadata["count_total_matches_metadata"].all()
        ),
        "source_counts_sha256": sha256_file(counts_path),
        "source_metadata_sha256": sha256_file(metadata_path),
        "note": (
            "GEO supplies sample/tissue labels but no spatial coordinates; this output supports "
            "patient-level differential analysis, not spot adjacency claims."
        ),
    }
    (results_dir / "GSE210037_spatial_pseudobulk_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == "__main__":
    main()
