#!/usr/bin/env python3
"""Build patient-level cohort and cell-type coverage tables for core irAE datasets."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


MIN_CELLS = 20
MIN_PATIENTS = 3


def clean_text(series: pd.Series) -> pd.Series:
    return series.astype("string").fillna("NA").replace({"nan": "NA", "<NA>": "NA"})


def unique_join(series: pd.Series) -> str:
    values = sorted({str(x) for x in series.dropna() if str(x) not in {"", "nan", "<NA>"}})
    return ";".join(values) if values else "NA"


def patient_inventory(
    obs: pd.DataFrame,
    dataset: str,
    organ: str,
    patient_col: str,
    condition_col: str,
    sample_col: str,
    celltype_col: str,
    extra_cols: list[str],
    analysis_role: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    keep = [patient_col, condition_col, sample_col, celltype_col] + [x for x in extra_cols if x in obs]
    frame = obs.loc[:, list(dict.fromkeys(keep))].copy()
    for col in frame.columns:
        frame[col] = clean_text(frame[col])

    patient_rows = []
    for (patient, condition), block in frame.groupby([patient_col, condition_col], observed=True, dropna=False):
        row = {
            "dataset": dataset,
            "organ": organ,
            "analysis_role": analysis_role,
            "patient_id": patient,
            "condition": condition,
            "n_cells": len(block),
            "n_samples": block[sample_col].nunique(dropna=True),
            "sample_ids": unique_join(block[sample_col]),
            "n_celltypes": block[celltype_col].nunique(dropna=True),
        }
        for col in extra_cols:
            if col in block:
                row[col] = unique_join(block[col])
        patient_rows.append(row)

    counts = (
        frame.groupby([patient_col, condition_col, celltype_col], observed=True, dropna=False)
        .size()
        .rename("n_cells")
        .reset_index()
        .rename(columns={patient_col: "patient_id", condition_col: "condition", celltype_col: "celltype"})
    )
    counts.insert(0, "organ", organ)
    counts.insert(0, "dataset", dataset)
    counts.insert(2, "analysis_role", analysis_role)
    counts["eligible_min_20_cells"] = counts["n_cells"] >= MIN_CELLS
    return pd.DataFrame(patient_rows), counts


def coverage_summary(counts: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (dataset, organ, analysis_role, condition, celltype), block in counts.groupby(
        ["dataset", "organ", "analysis_role", "condition", "celltype"], observed=True
    ):
        eligible = block.loc[block["eligible_min_20_cells"]]
        rows.append(
            {
                "dataset": dataset,
                "organ": organ,
                "analysis_role": analysis_role,
                "condition": condition,
                "celltype": celltype,
                "n_patients_total": block["patient_id"].nunique(),
                "n_patients_min_20_cells": eligible["patient_id"].nunique(),
                "median_cells_per_patient": float(block["n_cells"].median()),
                "min_cells_per_patient": int(block["n_cells"].min()),
                "max_cells_per_patient": int(block["n_cells"].max()),
                "passes_min_3_patients": eligible["patient_id"].nunique() >= MIN_PATIENTS,
            }
        )
    return pd.DataFrame(rows).sort_values(["dataset", "celltype", "condition"]).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True, type=Path)
    args = parser.parse_args()
    project = args.project.resolve()

    results_dir = project / "02_results"
    tables_dir = project / "tables_for_article"
    results_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    heart = pd.read_parquet(results_dir / "GSE228597_obs.parquet")
    colon = pd.read_parquet(results_dir / "GSE253720_obs.parquet")
    lung = pd.read_parquet(results_dir / "GSE277136_obs.parquet")
    lung["condition"] = lung["sample"].str.contains("-AE", regex=False).map(
        {True: "CIP", False: "healthy_control_external_source"}
    )

    heart_inventory, heart_counts = patient_inventory(
        heart,
        dataset="GSE228597",
        organ="heart",
        patient_col="donor",
        condition_col="condition",
        sample_col="Channel",
        celltype_col="lineage_names",
        extra_cols=[
            "fatal",
            "on_steroids",
            "combination_therapy",
            "gender",
            "donation_type",
            "region",
            "age_range",
            "institution",
        ],
        analysis_role="primary_with_institution_adjustment",
    )
    colon_inventory, colon_counts = patient_inventory(
        colon,
        dataset="GSE253720",
        organ="colon",
        patient_col="CoLabs_patient",
        condition_col="disease",
        sample_col="SAMPLE.by.SNPs",
        celltype_col="220414 JH COARSE",
        extra_cols=["cancer", "cpi", "cpi_sup", "suppression", "LIBRARY"],
        analysis_role="primary_discovery",
    )
    lung_inventory, lung_counts = patient_inventory(
        lung,
        dataset="GSE277136",
        organ="lung_BALF",
        patient_col="sample",
        condition_col="condition",
        sample_col="sample",
        celltype_col="cell_label",
        extra_cols=["source", "datasouce", "sub_label"],
        analysis_role="external_source_confounded_validation_only",
    )

    inventory = pd.concat([heart_inventory, colon_inventory, lung_inventory], ignore_index=True, sort=False)
    counts = pd.concat([heart_counts, colon_counts, lung_counts], ignore_index=True)
    coverage = coverage_summary(counts)

    inventory = inventory.sort_values(["dataset", "condition", "patient_id"]).reset_index(drop=True)
    counts = counts.sort_values(["dataset", "patient_id", "celltype"]).reset_index(drop=True)

    inventory.to_csv(tables_dir / "Table_S3_core_patient_inventory.tsv", sep="\t", index=False)
    coverage.to_csv(tables_dir / "Table_S4_celltype_patient_coverage.tsv", sep="\t", index=False)
    counts.to_csv(results_dir / "core_patient_celltype_counts.tsv.gz", sep="\t", index=False, compression="gzip")

    cohort = (
        inventory.groupby(["dataset", "organ", "analysis_role", "condition"], observed=True)
        .agg(n_patients=("patient_id", "nunique"), n_cells=("n_cells", "sum"), n_patient_condition_rows=("patient_id", "size"))
        .reset_index()
    )
    cohort.to_csv(results_dir / "core_cohort_summary.tsv", sep="\t", index=False)
    print(cohort.to_string(index=False))

    conflicting = inventory.groupby(["dataset", "patient_id"], observed=True)["condition"].nunique()
    conflicting = conflicting[conflicting > 1]
    print(f"\nPatients mapped to more than one condition: {len(conflicting)}")
    if len(conflicting):
        print(conflicting.to_string())

    passing = coverage.loc[coverage["passes_min_3_patients"]]
    print("\nCell-type/condition strata passing >=3 patients with >=20 cells:")
    print(passing.groupby(["dataset", "condition"], observed=True).size().rename("n_celltypes").to_string())


if __name__ == "__main__":
    main()
