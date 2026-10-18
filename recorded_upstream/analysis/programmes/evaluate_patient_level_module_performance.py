#!/usr/bin/env python
"""Evaluate frozen epithelial modules with patient/sample as the unit."""

from __future__ import annotations

import gzip
import pathlib
import sys
import zlib

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread


MODULES = [
    "ICI_COLITIS_EPITHELIAL_UP_STRICT",
    "ICI_COLITIS_EPITHELIAL_UP_NON_IFN",
    "ICI_COLITIS_EPITHELIAL_LOSS_STRICT",
    "CURATED_IFN_VISIBILITY",
]


def bh_adjust(pvalues: pd.Series) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    adj = np.minimum.accumulate((ranked * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
    restored = np.empty(len(pv))
    restored[order] = np.minimum(adj, 1.0)
    out[valid] = restored
    return out


def read_pb(root: pathlib.Path, prefix: str) -> tuple[sparse.csc_matrix, pd.DataFrame, pd.Index]:
    result = root / "02_results"
    meta = pd.read_csv(result / f"{prefix}_metadata.tsv", sep="\t")
    genes = pd.read_csv(result / f"{prefix}_genes.tsv", sep="\t")["gene"].astype(str)
    with gzip.open(result / f"{prefix}_counts.mtx.gz", "rb") as handle:
        matrix = mmread(handle)
    if matrix.shape == (len(meta), len(genes)):
        matrix = matrix.T
    if matrix.shape != (len(genes), len(meta)):
        raise RuntimeError(f"Unexpected {prefix} matrix shape {matrix.shape}")
    return matrix.tocsc(), meta, pd.Index(genes)


def aggregate_independent_families(
    matrix: sparse.csc_matrix, meta: pd.DataFrame, mapping: pd.DataFrame
) -> tuple[sparse.csc_matrix, pd.DataFrame]:
    family_map = dict(zip(mapping["cluster"].astype(str), mapping["cell_family"]))
    m = meta.copy()
    m["cluster_key"] = m["celltype"].astype(float).astype(int).astype(str)
    m["family"] = m["cluster_key"].map(family_map)
    keep = (
        m["condition"].isin(["irColitis", "On ICI therapy"])
        & m["drug"].eq("PD-1")
        & m["family"].isin(
            [
                "Absorptive epithelial",
                "Immature epithelial",
                "Mature absorptive epithelial",
                "Secretory epithelial",
            ]
        )
    )
    m = m.loc[keep].copy()
    original_indices = np.flatnonzero(keep.to_numpy())
    m["group"] = m["patient_id"].astype(str) + "||" + m["condition"] + "||" + m["family"]
    groups = pd.Index(m["group"].drop_duplicates())
    row = np.arange(len(m))
    col = groups.get_indexer(m["group"])
    membership = sparse.csc_matrix((np.ones(len(m)), (row, col)), shape=(len(m), len(groups)))
    aggregated = matrix[:, original_indices] @ membership
    records = []
    for group in groups:
        block = m.loc[m["group"].eq(group)]
        records.append(
            {
                "sample_id": group,
                "patient_id": block["patient_id"].iloc[0],
                "condition": block["condition"].iloc[0],
                "celltype": block["family"].iloc[0],
                "n_cells": int(block["n_cells"].sum()),
            }
        )
    return aggregated.tocsc(), pd.DataFrame(records)


def rank_scores(
    matrix: sparse.csc_matrix,
    genes: pd.Index,
    meta: pd.DataFrame,
    module_def: pd.DataFrame,
    dataset: str,
    cohort_role: str,
) -> pd.DataFrame:
    dense = matrix.toarray().astype(float)
    library = dense.sum(axis=0)
    if np.any(library <= 0):
        raise RuntimeError(f"Zero library in {dataset}")
    logcpm = np.log2(dense / library[None, :] * 1_000_000 + 1)
    ranks = np.apply_along_axis(stats.rankdata, 0, logcpm) / len(genes)
    symbol_to_index = {str(g).split("|")[-1]: i for i, g in enumerate(genes)}
    rows = []
    for module in MODULES:
        requested = module_def.loc[module_def["module"].eq(module), "gene"].astype(str).unique()
        indices = [symbol_to_index[g] for g in requested if g in symbol_to_index]
        if len(indices) < 10:
            raise RuntimeError(f"Only {len(indices)} genes available for {dataset}/{module}")
        raw = ranks[indices, :].mean(axis=0)
        expected = module_def.loc[module_def["module"].eq(module), "expected_direction_in_irAE"].iloc[0]
        oriented = raw if expected == "up" else 1 - raw
        for j, record in meta.reset_index(drop=True).iterrows():
            rows.append(
                {
                    "dataset": dataset,
                    "cohort_role": cohort_role,
                    "sample_id": record["sample_id"],
                    "patient_id": record["patient_id"],
                    "condition": record["condition"],
                    "celltype": record["celltype"],
                    "n_cells_or_spots": record["n_cells"],
                    "module": module,
                    "expected_direction_in_irAE": expected,
                    "n_genes_available": len(indices),
                    "raw_mean_within_sample_gene_rank": float(raw[j]),
                    "disease_oriented_score": float(oriented[j]),
                }
            )
    return pd.DataFrame(rows)


def bootstrap_auc(x: np.ndarray, y: np.ndarray, seed: int, n_boot: int = 5000) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values = np.empty(n_boot)
    for i in range(n_boot):
        xb = rng.choice(x, len(x), replace=True)
        yb = rng.choice(y, len(y), replace=True)
        values[i] = stats.mannwhitneyu(xb, yb, alternative="two-sided").statistic / (len(xb) * len(yb))
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def performance(scores: pd.DataFrame) -> pd.DataFrame:
    comparisons = {
        "GSE253720": ("CPI_colitis", "HC"),
        "GSE206300": ("irColitis", "On ICI therapy"),
        "GSE210037": ("irAE_colitis", "healthy_colon"),
    }
    rows = []
    for (dataset, role, celltype, module), block in scores.groupby(
        ["dataset", "cohort_role", "celltype", "module"], observed=True
    ):
        case, control = comparisons[dataset]
        x = block.loc[block["condition"].eq(case), "disease_oriented_score"].to_numpy()
        y = block.loc[block["condition"].eq(control), "disease_oriented_score"].to_numpy()
        if len(x) < 3 or len(y) < 3:
            continue
        test = stats.mannwhitneyu(x, y, alternative="two-sided")
        auc = float(test.statistic / (len(x) * len(y)))
        ci_low, ci_high = bootstrap_auc(
            x, y, seed=zlib.crc32(f"{dataset}|{celltype}|{module}".encode("utf-8"))
        )
        delta = float((np.greater.outer(x, y).sum() - np.less.outer(x, y).sum()) / (len(x) * len(y)))
        rows.append(
            {
                "dataset": dataset,
                "cohort_role": role,
                "celltype": celltype,
                "module": module,
                "case": case,
                "control": control,
                "n_case": len(x),
                "n_control": len(y),
                "median_case_oriented_score": float(np.median(x)),
                "median_control_oriented_score": float(np.median(y)),
                "cliffs_delta": delta,
                "auc_disease_higher": auc,
                "bootstrap_auc_ci_low": ci_low,
                "bootstrap_auc_ci_high": ci_high,
                "mannwhitney_p": float(test.pvalue),
            }
        )
    result = pd.DataFrame(rows)
    result["fdr_across_module_celltype_tests"] = bh_adjust(result["mannwhitney_p"])
    return result.sort_values(["cohort_role", "dataset", "celltype", "module"])


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: evaluate_patient_level_module_performance.py <project_dir>")
    root = pathlib.Path(sys.argv[1]).resolve()
    tables = root / "tables_for_article"
    module_def = pd.read_csv(tables / "Table_S21_frozen_module_definitions.tsv", sep="\t")
    all_scores = []

    matrix, meta, genes = read_pb(root, "GSE253720_pseudobulk")
    keep = meta["celltype"].eq("Epithelial") & meta["condition"].isin(["CPI_colitis", "HC"])
    dmeta = meta.loc[keep].copy()
    dmeta["sample_id"] = dmeta["pseudobulk_id"]
    dmeta = dmeta.rename(columns={"n_cells": "n_cells"})
    all_scores.append(
        rank_scores(
            matrix[:, np.flatnonzero(keep.to_numpy())], genes,
            dmeta[["sample_id", "patient_id", "condition", "celltype", "n_cells"]],
            module_def, "GSE253720", "module_derivation_self_check"
        )
    )

    matrix, meta, genes = read_pb(root, "GSE206300_pseudobulk")
    mapping = pd.read_csv(tables / "Table_S13_GSE206300_cluster_family_mapping.tsv", sep="\t")
    fam_matrix, fam_meta = aggregate_independent_families(matrix, meta, mapping)
    all_scores.append(
        rank_scores(
            fam_matrix, genes, fam_meta, module_def,
            "GSE206300", "independent_primary_PD1_only"
        )
    )

    matrix, meta, genes = read_pb(root, "GSE210037_spatial_pseudobulk")
    keep = meta["tissue"].isin(["irAE_colitis", "healthy_colon"])
    smeta = meta.loc[keep].copy()
    smeta["sample_id"] = smeta["sample"]
    smeta["patient_id"] = smeta["sample"]
    smeta["condition"] = smeta["tissue"]
    smeta["celltype"] = "whole_spot_sample"
    smeta["n_cells"] = smeta["n_spots"]
    all_scores.append(
        rank_scores(
            matrix[:, np.flatnonzero(keep.to_numpy())], genes,
            smeta[["sample_id", "patient_id", "condition", "celltype", "n_cells"]],
            module_def, "GSE210037", "independent_spatial_platform"
        )
    )

    scores = pd.concat(all_scores, ignore_index=True)
    perf = performance(scores)
    scores.to_csv(tables / "Table_S55_patient_level_frozen_module_scores.tsv", sep="\t", index=False)
    perf.to_csv(
        tables / "Table_S56_patient_level_frozen_module_classification_performance.tsv",
        sep="\t", index=False
    )
    print("Patient-level frozen-module performance complete")
    print(
        perf.loc[perf["cohort_role"].ne("module_derivation_self_check"),
                 ["dataset", "celltype", "module", "auc_disease_higher", "mannwhitney_p",
                  "fdr_across_module_celltype_tests"]].to_string(index=False)
    )


if __name__ == "__main__":
    main()
