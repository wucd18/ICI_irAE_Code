#!/usr/bin/env python
"""Patient-replicated localization of the TNFSF12--TNFRSF12A--SORL1 axis.

The analysis uses cell-level counts for descriptive localization, but all
case-control tests use patient-level summaries to avoid cell-level
pseudoreplication.
"""

from __future__ import annotations

import pathlib
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse, stats


GENES = ["TNFSF12", "TNFRSF12A", "SORL1", "DSC2", "PLPP1"]


def bh_adjust(values: pd.Series) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return out
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    adjusted = np.minimum.accumulate((ranked * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
    restored = np.empty(len(pv))
    restored[order] = np.minimum(adjusted, 1.0)
    out[valid] = restored
    return out


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    return float((np.greater.outer(x, y).sum() - np.less.outer(x, y).sum()) / (len(x) * len(y)))


def select_gene_indices(var_names: pd.Index) -> tuple[list[int], list[str]]:
    symbols = pd.Index([str(x).split("|")[-1] for x in var_names])
    indices: list[int] = []
    retained: list[str] = []
    for gene in GENES:
        hit = np.flatnonzero(symbols == gene)
        if len(hit) == 1:
            indices.append(int(hit[0]))
            retained.append(gene)
    if "SORL1" not in retained or "TNFRSF12A" not in retained:
        raise RuntimeError(f"Required axis genes not found; retained={retained}")
    return indices, retained


def patient_gene_summaries(
    matrix,
    obs: pd.DataFrame,
    dataset: str,
    patient_col: str,
    condition_col: str,
    celltype_col: str,
    gene_names: list[str],
    extra_cols: list[str] | None = None,
) -> pd.DataFrame:
    if sparse.issparse(matrix):
        matrix = matrix.tocsr()
    else:
        matrix = np.asarray(matrix)
    extra_cols = extra_cols or []
    rows: list[dict] = []
    group_cols = [patient_col, condition_col, celltype_col] + extra_cols
    for key, index in obs.groupby(group_cols, observed=True, sort=False).indices.items():
        if not isinstance(key, tuple):
            key = (key,)
        idx = np.asarray(index, dtype=int)
        n_cells = len(idx)
        if n_cells < 20:
            continue
        block = matrix[idx, :]
        if sparse.issparse(block):
            detected = np.asarray((block > 0).sum(axis=0)).ravel()
            means = np.asarray(block.mean(axis=0)).ravel()
            mean_log1p = np.asarray(block.copy().astype(float).log1p().mean(axis=0)).ravel()
        else:
            detected = (block > 0).sum(axis=0)
            means = block.mean(axis=0)
            mean_log1p = np.log1p(block).mean(axis=0)
        base = dict(zip(group_cols, key))
        for j, gene in enumerate(gene_names):
            rows.append(
                {
                    "dataset": dataset,
                    "patient_id": str(base[patient_col]),
                    "condition": str(base[condition_col]),
                    "celltype": str(base[celltype_col]),
                    **{col: str(base[col]) for col in extra_cols},
                    "gene": gene,
                    "n_cells": n_cells,
                    "detected_cells": int(detected[j]),
                    "detected_fraction": float(detected[j] / n_cells),
                    "mean_raw_count_per_cell": float(means[j]),
                    "mean_log1p_raw_count_per_cell": float(mean_log1p[j]),
                }
            )
    return pd.DataFrame(rows)


def compare_conditions(
    patient: pd.DataFrame,
    case: str,
    control: str,
    analysis_name: str,
) -> pd.DataFrame:
    rows: list[dict] = []
    for (dataset, celltype, gene), block in patient.groupby(
        ["dataset", "celltype", "gene"], observed=True, sort=False
    ):
        x = block.loc[block["condition"] == case, "mean_log1p_raw_count_per_cell"].dropna().to_numpy()
        y = block.loc[block["condition"] == control, "mean_log1p_raw_count_per_cell"].dropna().to_numpy()
        if len(x) < 3 or len(y) < 3:
            continue
        test = stats.mannwhitneyu(x, y, alternative="two-sided")
        rows.append(
            {
                "analysis": analysis_name,
                "dataset": dataset,
                "celltype": celltype,
                "gene": gene,
                "case": case,
                "control": control,
                "n_case_patients": len(x),
                "n_control_patients": len(y),
                "median_case": float(np.median(x)),
                "median_control": float(np.median(y)),
                "median_difference": float(np.median(x) - np.median(y)),
                "cliffs_delta": cliffs_delta(x, y),
                "p": float(test.pvalue),
            }
        )
    result = pd.DataFrame(rows)
    if len(result):
        result["fdr"] = bh_adjust(result["p"])
    return result


def summarize_localization(patient: pd.DataFrame) -> pd.DataFrame:
    return (
        patient.groupby(["dataset", "condition", "celltype", "gene"], observed=True)
        .agg(
            n_patients=("patient_id", "nunique"),
            total_cells=("n_cells", "sum"),
            median_patient_detected_fraction=("detected_fraction", "median"),
            q1_patient_detected_fraction=("detected_fraction", lambda x: x.quantile(0.25)),
            q3_patient_detected_fraction=("detected_fraction", lambda x: x.quantile(0.75)),
            median_patient_mean_log1p_count=("mean_log1p_raw_count_per_cell", "median"),
        )
        .reset_index()
    )


def discovery(root: pathlib.Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = root / "01_data" / "GSE253720_Biopsy_RNA_Final.h5ad"
    a = ad.read_h5ad(path, backed="r")
    try:
        indices, retained = select_gene_indices(a.var_names)
        obs = a.obs[["CoLabs_patient", "disease", "220414 JH COARSE"]].copy().reset_index(drop=True)
        obs.columns = ["patient_id", "condition", "celltype"]
        matrix = a.layers["counts"][:, indices]
        patient = patient_gene_summaries(
            matrix, obs, "GSE253720", "patient_id", "condition", "celltype", retained
        )
    finally:
        a.file.close()
    comparison = compare_conditions(
        patient, "CPI_colitis", "HC", "CPI_colitis_vs_HC_patient_level"
    )
    return patient, comparison


def independent(root: pathlib.Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = root / "01_data" / "GSE206300_ircolitis-tissue-epithelial.h5ad"
    mapping = pd.read_csv(
        root / "tables_for_article" / "Table_S13_GSE206300_cluster_family_mapping.tsv",
        sep="\t",
        dtype={"cluster": str},
    )
    family_map = dict(zip(mapping["cluster"], mapping["cell_family"]))
    a = ad.read_h5ad(path, backed="r")
    try:
        indices, retained = select_gene_indices(a.var_names)
        obs = a.obs[["donor", "class", "drug", "cluster"]].copy().reset_index(drop=True)
        obs["cluster"] = obs["cluster"].astype(float).astype(int).astype(str)
        obs["cell_family"] = obs["cluster"].map(family_map)
        obs = obs.loc[obs["cell_family"].notna()].copy()
        # The primary independent contrast is deliberately restricted to PD-1.
        keep = obs["class"].isin(["irColitis", "On ICI therapy"]) & obs["drug"].eq("PD-1")
        retained_rows = np.flatnonzero(keep.to_numpy())
        selected_obs = obs.loc[keep, ["donor", "class", "cell_family", "drug"]].reset_index(drop=True)
        matrix = a.X[retained_rows, :][:, indices]
        patient = patient_gene_summaries(
            matrix,
            selected_obs,
            "GSE206300",
            "donor",
            "class",
            "cell_family",
            retained,
            extra_cols=["drug"],
        )
    finally:
        a.file.close()
    comparison = compare_conditions(
        patient,
        "irColitis",
        "On ICI therapy",
        "irColitis_vs_ICI_control_PD1_only_patient_level",
    )
    return patient, comparison


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: analyze_TWEAK_SORL1_axis.py <project_dir>")
    root = pathlib.Path(sys.argv[1]).resolve()
    tables = root / "tables_for_article"
    d_patient, d_comp = discovery(root)
    i_patient, i_comp = independent(root)
    patient = pd.concat([d_patient, i_patient], ignore_index=True, sort=False)
    localization = summarize_localization(patient)
    comparisons = pd.concat([d_comp, i_comp], ignore_index=True, sort=False)
    patient.to_csv(tables / "Table_S44_TWEAK_SORL1_axis_patient_expression.tsv", sep="\t", index=False)
    localization.to_csv(
        tables / "Table_S45_TWEAK_SORL1_axis_celltype_localization.tsv", sep="\t", index=False
    )
    comparisons.to_csv(
        tables / "Table_S46_TWEAK_SORL1_axis_patient_level_comparisons.tsv", sep="\t", index=False
    )
    focus = comparisons.loc[comparisons["gene"].isin(["TNFRSF12A", "SORL1"])].copy()
    print("TWEAK--Fn14--SORL1 axis analysis complete")
    print(focus.sort_values(["dataset", "gene", "p"]).to_string(index=False))


if __name__ == "__main__":
    main()
