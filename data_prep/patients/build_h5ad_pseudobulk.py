#!/usr/bin/env python3
"""Aggregate raw single-cell counts to patient-level pseudobulk matrices."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


def first_nonmissing(values: pd.Series) -> str:
    vals = sorted({str(x) for x in values.dropna() if str(x) not in {"", "nan", "<NA>"}})
    return ";".join(vals) if vals else "NA"


def get_matrix_chunk(adata: ad.AnnData, source: str, start: int, stop: int):
    if source == "raw":
        if adata.raw is None:
            raise ValueError("Requested raw counts, but adata.raw is absent")
        return adata.raw.X[start:stop]
    if source.startswith("layer:"):
        layer = source.split(":", 1)[1]
        if layer not in adata.layers:
            raise ValueError(f"Layer {layer!r} is absent")
        return adata.layers[layer][start:stop]
    if source == "X":
        return adata.X[start:stop]
    if source == "raw_log1p_1e4":
        if adata.raw is None:
            raise ValueError("Requested reconstructed raw counts, but adata.raw is absent")
        if "total_counts" not in adata.obs:
            raise ValueError("raw_log1p_1e4 requires obs['total_counts']")
        matrix = sparse.csr_matrix(adata.raw.X[start:stop]).astype(np.float64)
        matrix.data = np.expm1(matrix.data)
        totals = adata.obs.iloc[start:stop]["total_counts"].to_numpy(dtype=np.float64)
        matrix = sparse.diags(totals / 10000.0) @ matrix
        deviation = np.max(np.abs(matrix.data - np.rint(matrix.data))) if matrix.nnz else 0.0
        if deviation > 0.01:
            raise ValueError(f"Log-normalized count reconstruction failed; max deviation={deviation}")
        matrix.data = np.rint(matrix.data)
        matrix.eliminate_zeros()
        return matrix
    raise ValueError(f"Unsupported count source: {source}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output_prefix", type=Path)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--organ", required=True)
    parser.add_argument("--patient-col", required=True)
    parser.add_argument("--condition-col")
    parser.add_argument("--condition-from-sample", action="store_true")
    parser.add_argument("--celltype-col", required=True)
    parser.add_argument(
        "--count-source",
        required=True,
        help="raw, X, layer:<name>, or raw_log1p_1e4 (exact reconstruction using total_counts)",
    )
    parser.add_argument("--strata-cols", default="", help="Comma-separated columns that must not be merged")
    parser.add_argument("--extra-cols", default="", help="Comma-separated metadata columns to summarize")
    parser.add_argument("--chunk-size", type=int, default=2000)
    args = parser.parse_args()

    out = args.output_prefix.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.input.resolve(), backed="r")
    obs = adata.obs.copy()

    if args.condition_from_sample:
        obs["_derived_condition"] = np.where(
            obs[args.patient_col].astype(str).str.contains("AE", case=False, regex=False),
            "CIP",
            "control",
        )
        condition_col = "_derived_condition"
    else:
        if not args.condition_col:
            raise ValueError("--condition-col is required unless --condition-from-sample is used")
        condition_col = args.condition_col

    strata = [x for x in args.strata_cols.split(",") if x]
    extras = [x for x in args.extra_cols.split(",") if x]
    required = [args.patient_col, condition_col, args.celltype_col] + strata
    missing = [x for x in required if x not in obs.columns]
    if missing:
        raise KeyError(f"Missing required obs columns: {missing}")

    group_cols = [args.patient_col, condition_col, args.celltype_col] + strata
    key_frame = obs[group_cols].astype("string").fillna("NA")
    key_frame.columns = ["patient_id", "condition", "celltype"] + strata
    group_index = pd.MultiIndex.from_frame(key_frame)
    codes, uniques = pd.factorize(group_index, sort=True)
    n_groups = len(uniques)

    if args.count_source in {"raw", "raw_log1p_1e4"}:
        genes = pd.Index(adata.raw.var_names.astype(str), name="gene")
    else:
        genes = pd.Index(adata.var_names.astype(str), name="gene")

    aggregate = sparse.csr_matrix((n_groups, len(genes)), dtype=np.float64)
    n_cells = np.bincount(codes, minlength=n_groups)
    for start in range(0, adata.n_obs, args.chunk_size):
        stop = min(start + args.chunk_size, adata.n_obs)
        matrix = get_matrix_chunk(adata, args.count_source, start, stop)
        matrix = sparse.csr_matrix(matrix)
        membership = sparse.csr_matrix(
            (np.ones(stop - start), (codes[start:stop], np.arange(stop - start))),
            shape=(n_groups, stop - start),
        )
        aggregate = aggregate + membership @ matrix
        print(f"{args.dataset}: {stop}/{adata.n_obs} cells", flush=True)

    if aggregate.nnz:
        max_fractional = float(np.max(np.abs(aggregate.data - np.rint(aggregate.data))))
    else:
        max_fractional = 0.0
    if max_fractional > 1e-6:
        raise ValueError(f"Aggregated values are not integer counts; max fractional deviation={max_fractional}")
    aggregate.data = np.rint(aggregate.data).astype(np.int64)
    aggregate.eliminate_zeros()

    meta = uniques.to_frame(index=False)
    meta.columns = ["patient_id", "condition", "celltype"] + strata
    meta.insert(0, "organ", args.organ)
    meta.insert(0, "dataset", args.dataset)
    meta["n_cells"] = n_cells
    meta["eligible_min_20_cells"] = meta["n_cells"] >= 20
    meta.insert(0, "pseudobulk_id", [f"{args.dataset}_PB{i + 1:04d}" for i in range(n_groups)])

    obs_with_code = obs.copy()
    obs_with_code["_group_code"] = codes
    for col in extras:
        if col in obs_with_code:
            summarized = obs_with_code.groupby("_group_code", observed=True)[col].agg(first_nonmissing)
            meta[col] = meta.index.to_series().map(summarized).fillna("NA")

    meta.to_csv(out.with_name(out.name + "_metadata.tsv"), sep="\t", index=False)
    pd.DataFrame({"gene_index": np.arange(len(genes)), "gene": genes}).to_csv(
        out.with_name(out.name + "_genes.tsv"), sep="\t", index=False
    )

    mtx_plain = out.with_name(out.name + "_counts.mtx")
    mtx_gz = out.with_name(out.name + "_counts.mtx.gz")
    with mtx_plain.open("wb") as handle:
        mmwrite(handle, aggregate, field="integer", symmetry="general")
    with mtx_plain.open("rb") as src, gzip.open(mtx_gz, "wb", compresslevel=6) as dst:
        shutil.copyfileobj(src, dst, length=8 * 1024 * 1024)
    mtx_plain.unlink()

    audit = {
        "dataset": args.dataset,
        "organ": args.organ,
        "input": str(args.input.resolve()),
        "n_cells": int(adata.n_obs),
        "n_genes": int(len(genes)),
        "n_pseudobulks": int(n_groups),
        "count_source": args.count_source,
        "group_columns": group_cols,
        "max_fractional_deviation": max_fractional,
        "matrix_nnz": int(aggregate.nnz),
    }
    out.with_name(out.name + "_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
