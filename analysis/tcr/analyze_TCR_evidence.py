#!/usr/bin/env python3
"""Patient-level TCR evidence for ICI colitis and myocarditis decoupling."""

from __future__ import annotations

import argparse
import gzip
import math
import pathlib
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread


SEED = 20260826


def bh_adjust(values: pd.Series) -> pd.Series:
    arr = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    out = np.full(len(arr), np.nan)
    valid = np.isfinite(arr)
    p = arr[valid]
    if not len(p):
        return pd.Series(out, index=values.index)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    out[valid] = restored
    return pd.Series(out, index=values.index)


def cliffs_delta(group1: np.ndarray, group0: np.ndarray) -> float:
    group1 = np.asarray(group1, dtype=float)
    group0 = np.asarray(group0, dtype=float)
    comparisons = group1[:, None] - group0[None, :]
    return float((np.sum(comparisons > 0) - np.sum(comparisons < 0)) / comparisons.size)


def repertoire_metrics(counts: pd.Series) -> dict[str, float]:
    counts = pd.to_numeric(counts, errors="coerce").dropna()
    counts = counts[counts > 0]
    total = float(counts.sum())
    richness = int(len(counts))
    if total <= 0 or richness == 0:
        return {}
    p = counts.to_numpy(dtype=float) / total
    shannon = float(-np.sum(p * np.log(p)))
    normalized_clonality = float(1 - shannon / np.log(richness)) if richness > 1 else 1.0
    sorted_p = np.sort(p)[::-1]
    return {
        "n_tcr_cells": int(total),
        "n_clonotypes": richness,
        "expanded_cell_fraction_ge2": float(counts[counts >= 2].sum() / total),
        "large_clone_cell_fraction_ge5": float(counts[counts >= 5].sum() / total),
        "shannon_entropy": shannon,
        "normalized_clonality": normalized_clonality,
        "simpson_concentration": float(np.sum(p**2)),
        "top1_fraction": float(sorted_p[:1].sum()),
        "top10_fraction": float(sorted_p[:10].sum()),
    }


def rarefy_metrics(counts: pd.Series, depth: int, iterations: int, seed: int) -> dict[str, float]:
    counts = pd.to_numeric(counts, errors="coerce").dropna().astype(int)
    counts = counts[counts > 0]
    total = int(counts.sum())
    if total < depth or depth < 2:
        return {
            f"rarefied_richness_{depth}": np.nan,
            f"rarefied_clonality_{depth}": np.nan,
            f"rarefied_top10_fraction_{depth}": np.nan,
        }
    labels = np.repeat(np.arange(len(counts)), counts.to_numpy())
    rng = np.random.default_rng(seed)
    richness_values: list[float] = []
    clonality_values: list[float] = []
    top10_values: list[float] = []
    for _ in range(iterations):
        sampled = rng.choice(labels, size=depth, replace=False)
        sample_counts = pd.Series(sampled).value_counts()
        metrics = repertoire_metrics(sample_counts)
        richness_values.append(metrics["n_clonotypes"])
        clonality_values.append(metrics["normalized_clonality"])
        top10_values.append(metrics["top10_fraction"])
    return {
        f"rarefied_richness_{depth}": float(np.mean(richness_values)),
        f"rarefied_clonality_{depth}": float(np.mean(clonality_values)),
        f"rarefied_top10_fraction_{depth}": float(np.mean(top10_values)),
    }


def patient_clone_metrics(
    cells: pd.DataFrame,
    patient_col: str,
    clone_col: str,
    group_cols: list[str],
    rarefaction_depth: int,
) -> pd.DataFrame:
    clone_counts = (
        cells.groupby(group_cols + [patient_col, clone_col], observed=True)
        .size()
        .rename("clone_cells")
        .reset_index()
    )
    rows: list[dict] = []
    grouping = group_cols + [patient_col]
    for keys, block in clone_counts.groupby(grouping, observed=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(grouping, keys))
        row.update(repertoire_metrics(block.set_index(clone_col)["clone_cells"]))
        row.update(
            rarefy_metrics(
                block.set_index(clone_col)["clone_cells"],
                depth=rarefaction_depth,
                iterations=500,
                seed=SEED + sum(ord(char) for char in str(keys)),
            )
        )
        rows.append(row)
    return pd.DataFrame(rows)


def group_tests(
    table: pd.DataFrame,
    group_col: str,
    group1: str,
    group0: str,
    metrics: list[str],
    analysis: str,
    minimum_per_group: int,
) -> pd.DataFrame:
    rows: list[dict] = []
    for metric in metrics:
        x1 = pd.to_numeric(table.loc[table[group_col] == group1, metric], errors="coerce").dropna().to_numpy()
        x0 = pd.to_numeric(table.loc[table[group_col] == group0, metric], errors="coerce").dropna().to_numpy()
        if len(x1) < minimum_per_group or len(x0) < minimum_per_group:
            continue
        test = stats.mannwhitneyu(x1, x0, alternative="two-sided", method="auto")
        rows.append(
            {
                "analysis": analysis,
                "metric": metric,
                "group1": group1,
                "group0": group0,
                "n_group1": len(x1),
                "n_group0": len(x0),
                "median_group1": float(np.median(x1)),
                "median_group0": float(np.median(x0)),
                "median_difference": float(np.median(x1) - np.median(x0)),
                "cliffs_delta": cliffs_delta(x1, x0),
                "p": float(test.pvalue),
            }
        )
    result = pd.DataFrame(rows)
    if len(result):
        result["fdr_within_analysis"] = bh_adjust(result["p"])
    return result


def analyze_colon_tcr(results_dir: pathlib.Path, tables_dir: pathlib.Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    obs = pd.read_parquet(results_dir / "GSE253720_obs.parquet")
    coarse = "220414 JH COARSE"
    fine = "220414 JH FINE"
    valid_lineage = obs[coarse].astype(str).isin(["CD4T", "CD8T"]) | (obs[fine].astype(str) == "Cycling T")
    valid_tcr = ~obs["TCR.cdr3s_aa"].astype(str).isin(["NA", "nan", "None", ""])
    cells = obs.loc[
        valid_lineage & valid_tcr,
        ["cell_barcode", "CoLabs_patient", "disease", coarse, fine, "TCR.cdr3s_aa"],
    ].copy()
    cells = cells.rename(
        columns={
            "CoLabs_patient": "patient_id",
            coarse: "coarse_celltype",
            fine: "fine_celltype",
            "TCR.cdr3s_aa": "clone_key",
        }
    )
    cells["coarse_celltype"] = cells["coarse_celltype"].astype(str)
    cells["fine_celltype"] = cells["fine_celltype"].astype(str)
    cells["clone_key"] = cells["clone_key"].astype(str)

    patient_metrics = patient_clone_metrics(
        cells, "patient_id", "clone_key", ["disease"], rarefaction_depth=200
    )

    clone_states = (
        cells.groupby(["patient_id", "clone_key"], observed=True)
        .agg(
            disease=("disease", "first"),
            clone_cells=("cell_barcode", "size"),
            coarse_states=("coarse_celltype", lambda x: set(x)),
            fine_states=("fine_celltype", lambda x: set(x)),
        )
        .reset_index()
    )
    clone_states["shared_cd4_cd8"] = clone_states["coarse_states"].map(
        lambda states: {"CD4T", "CD8T"}.issubset(states)
    )
    cd4_resident = {"CD4T RM", "CD4T memory"}
    cd8_effector = {"CD8T RM", "CD8T GZMB", "CD8T GZMK", "Cycling T"}
    clone_states["shared_resident_effector"] = clone_states["fine_states"].map(
        lambda states: bool(states & cd4_resident) and bool(states & cd8_effector)
    )
    sharing = clone_states.groupby("patient_id", observed=True).apply(
        lambda block: pd.Series(
            {
                "cross_cd4_cd8_clone_fraction": float(block["shared_cd4_cd8"].mean()),
                "cross_cd4_cd8_cell_fraction": float(
                    block.loc[block["shared_cd4_cd8"], "clone_cells"].sum() / block["clone_cells"].sum()
                ),
                "resident_effector_shared_clone_fraction": float(block["shared_resident_effector"].mean()),
                "resident_effector_shared_cell_fraction": float(
                    block.loc[block["shared_resident_effector"], "clone_cells"].sum()
                    / block["clone_cells"].sum()
                ),
            }
        ),
        include_groups=False,
    ).reset_index()
    patient_metrics = patient_metrics.merge(sharing, on="patient_id", how="left")
    patient_metrics.to_csv(tables_dir / "Table_S31_colon_TCR_patient_metrics.tsv", sep="\t", index=False)

    clone_size = cells.groupby(["patient_id", "clone_key"], observed=True).size().rename("clone_cells")
    cells = cells.join(clone_size, on=["patient_id", "clone_key"])
    cells["expanded"] = cells["clone_cells"] >= 2
    state = (
        cells.groupby(["disease", "patient_id", "fine_celltype"], observed=True)
        .agg(
            n_tcr_cells=("cell_barcode", "size"),
            n_expanded_cells=("expanded", "sum"),
            n_clonotypes=("clone_key", "nunique"),
        )
        .reset_index()
    )
    state["expanded_cell_fraction"] = state["n_expanded_cells"] / state["n_tcr_cells"]
    state.to_csv(tables_dir / "Table_S32_colon_TCR_state_expansion.tsv", sep="\t", index=False)

    overall_metrics = [
        "expanded_cell_fraction_ge2",
        "large_clone_cell_fraction_ge5",
        "normalized_clonality",
        "simpson_concentration",
        "top10_fraction",
        "rarefied_clonality_200",
        "rarefied_top10_fraction_200",
        "cross_cd4_cd8_cell_fraction",
        "resident_effector_shared_cell_fraction",
    ]
    overall_tests = group_tests(
        patient_metrics,
        "disease",
        "CPI_colitis",
        "HC",
        overall_metrics,
        "colon_CPI_colitis_vs_HC_patient_level",
        minimum_per_group=3,
    )

    state_tests: list[pd.DataFrame] = []
    for state_name, block in state.groupby("fine_celltype", observed=True):
        eligible = block[block["n_tcr_cells"] >= 20]
        test = group_tests(
            eligible,
            "disease",
            "CPI_colitis",
            "HC",
            ["expanded_cell_fraction"],
            f"colon_state_{state_name}_CPI_colitis_vs_HC",
            minimum_per_group=3,
        )
        if len(test):
            test["cell_state"] = state_name
            state_tests.append(test)
    state_test_table = pd.concat(state_tests, ignore_index=True) if state_tests else pd.DataFrame()
    if len(state_test_table):
        state_test_table["fdr_across_states"] = bh_adjust(state_test_table["p"])
    return overall_tests, state_test_table


def couple_colon_tcr_to_epithelial_modules(
    root: pathlib.Path, results_dir: pathlib.Path, tables_dir: pathlib.Path
) -> pd.DataFrame:
    patient_metrics = pd.read_csv(tables_dir / "Table_S31_colon_TCR_patient_metrics.tsv", sep="\t")
    pb_meta = pd.read_csv(results_dir / "GSE253720_pseudobulk_metadata.tsv", sep="\t")
    genes = pd.read_csv(results_dir / "GSE253720_pseudobulk_genes.tsv", sep="\t")["gene"].astype(str)
    with gzip.open(results_dir / "GSE253720_pseudobulk_counts.mtx.gz", "rb") as handle:
        counts = mmread(handle)
    if counts.shape == (len(pb_meta), len(genes)):
        counts = counts.T
    elif counts.shape != (len(genes), len(pb_meta)):
        raise RuntimeError(f"Unexpected GSE253720 pseudobulk shape: {counts.shape}")
    counts = counts.tocsc()
    epithelial_indices = np.flatnonzero(
        (pb_meta["celltype"].astype(str) == "Epithelial")
        & pb_meta["eligible_min_20_cells"].astype(str).str.lower().isin(["true", "1"])
    )
    epithelial_meta = pb_meta.iloc[epithelial_indices].reset_index(drop=True)
    epithelial_counts = counts[:, epithelial_indices].toarray().astype(float)
    library_size = epithelial_counts.sum(axis=0)
    logcpm = np.log2(epithelial_counts / library_size[None, :] * 1_000_000 + 1)
    ranks = np.apply_along_axis(stats.rankdata, 0, logcpm) / logcpm.shape[0]

    module_def = pd.read_csv(tables_dir / "Table_S21_frozen_module_definitions.tsv", sep="\t")
    selected_modules = [
        "ICI_COLITIS_EPITHELIAL_UP_STRICT",
        "ICI_COLITIS_EPITHELIAL_UP_NON_IFN",
        "ICI_COLITIS_EPITHELIAL_LOSS_STRICT",
        "CURATED_IFN_VISIBILITY",
    ]
    gene_to_index = {gene: index for index, gene in enumerate(genes)}
    score_table = epithelial_meta[["patient_id", "condition", "n_cells"]].copy()
    for module in selected_modules:
        module_genes = module_def.loc[module_def["module"] == module, "gene"].astype(str)
        available = [gene_to_index[gene] for gene in module_genes if gene in gene_to_index]
        if len(available) < 10:
            raise RuntimeError(f"Insufficient coverage for {module}: {len(available)} genes")
        score_table[module] = ranks[available, :].mean(axis=0)
        score_table[f"{module}_n_genes"] = len(available)
    coupled = patient_metrics.merge(score_table, on="patient_id", how="inner", validate="one_to_one")
    coupled.to_csv(
        tables_dir / "Table_S36_colon_TCR_epithelial_module_patient_coupling.tsv",
        sep="\t",
        index=False,
    )

    tcr_metrics = [
        "expanded_cell_fraction_ge2",
        "normalized_clonality",
        "cross_cd4_cd8_cell_fraction",
        "resident_effector_shared_cell_fraction",
    ]
    cpi = coupled.loc[coupled["disease"] == "CPI_colitis"].copy()
    rows: list[dict] = []
    for tcr_metric in tcr_metrics:
        for module in selected_modules:
            valid = cpi[[tcr_metric, module]].dropna()
            if len(valid) < 7:
                continue
            result = stats.spearmanr(valid[tcr_metric], valid[module], alternative="two-sided")
            rows.append(
                {
                    "analysis": "within_CPI_colitis_patient_exploratory",
                    "n_patients": len(valid),
                    "tcr_metric": tcr_metric,
                    "epithelial_module": module,
                    "spearman_rho": float(result.statistic),
                    "p": float(result.pvalue),
                }
            )
    correlations = pd.DataFrame(rows)
    correlations["fdr"] = bh_adjust(correlations["p"])
    correlations.to_csv(
        tables_dir / "Table_S37_colon_TCR_epithelial_module_correlations.tsv",
        sep="\t",
        index=False,
    )
    return correlations


def cardiac_cell_clone_table(data_dir: pathlib.Path, results_dir: pathlib.Path) -> pd.DataFrame:
    obs = pd.read_parquet(results_dir / "GSE228597_obs.parquet")
    obs = obs.loc[
        obs["lineage_names"].astype(str) == "T and NK cells",
        ["cell_barcode", "Channel", "donor", "condition"],
    ].copy()
    obs["cell_barcode"] = obs["cell_barcode"].astype(str)
    cell_rows: list[dict] = []
    pattern = re.compile(r"^GSM\d+_(.+)_all_contig_annotations_tcr\.csv\.gz$")
    for path in sorted(data_dir.glob("GSM7133*_heart_all_contig_annotations_tcr.csv.gz")):
        match = pattern.match(path.name)
        if not match:
            continue
        channel = match.group(1)
        contigs = pd.read_csv(path)
        productive = contigs["productive"].fillna(False).astype(bool)
        confident = contigs["high_confidence"].fillna(False).astype(bool)
        valid = productive & confident & contigs["chain"].isin(["TRA", "TRB"]) & contigs["cdr3"].notna()
        contigs = contigs.loc[valid, ["barcode", "chain", "cdr3"]].copy()
        for barcode, block in contigs.groupby("barcode", observed=True):
            trb = sorted(set(block.loc[block["chain"] == "TRB", "cdr3"].astype(str)))
            if not trb:
                continue
            tra = sorted(set(block.loc[block["chain"] == "TRA", "cdr3"].astype(str)))
            clone_key = "TRA:" + ",".join(tra) + "|TRB:" + ",".join(trb)
            clean_barcode = str(barcode).removesuffix("-1")
            cell_rows.append(
                {
                    "cell_barcode": f"{channel}_{clean_barcode}",
                    "Channel": channel,
                    "clone_key": clone_key,
                    "trb_cdr3": ";".join(trb),
                }
            )
    cells = pd.DataFrame(cell_rows)
    cells = cells.merge(obs, on=["cell_barcode", "Channel"], how="inner", validate="one_to_one")
    return cells


def parse_soft_samples(path: pathlib.Path, gsm_prefix: str) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.read().splitlines()
    records: list[dict] = []
    current: dict | None = None
    for line in lines:
        if line.startswith("^SAMPLE = "):
            if current and current["gsm"].startswith(gsm_prefix):
                records.append(current)
            current = {"gsm": line.split("=", 1)[1].strip()}
        elif current is not None and line.startswith("!Sample_title = "):
            current["title"] = line.split("=", 1)[1].strip()
        elif current is not None and line.startswith("!Sample_source_name_ch1 = "):
            current["source_name"] = line.split("=", 1)[1].strip()
        elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
            value = line.split("=", 1)[1].strip()
            if ":" in value:
                key, item = value.split(":", 1)
                current[re.sub(r"[^a-z0-9]+", "_", key.lower()).strip("_")] = item.strip()
    if current and current["gsm"].startswith(gsm_prefix):
        records.append(current)
    return pd.DataFrame(records)


def morisita_horn(a: pd.Series, b: pd.Series) -> float:
    joined = pd.concat([a.rename("a"), b.rename("b")], axis=1).fillna(0.0)
    pa = joined["a"].to_numpy(dtype=float).copy()
    pb = joined["b"].to_numpy(dtype=float).copy()
    if pa.sum() <= 0 or pb.sum() <= 0:
        return np.nan
    pa /= pa.sum()
    pb /= pb.sum()
    denominator = np.sum(pa**2) + np.sum(pb**2)
    return float(2 * np.sum(pa * pb) / denominator) if denominator > 0 else np.nan


def overlap_metrics(a: pd.Series, b: pd.Series, label_a: str, label_b: str) -> dict[str, float]:
    a = a[a > 0]
    b = b[b > 0]
    shared = a.index.intersection(b.index)
    union = a.index.union(b.index)
    top_a = set(a.sort_values(ascending=False).head(20).index)
    top_b = set(b.sort_values(ascending=False).head(20).index)
    return {
        "compartment_a": label_a,
        "compartment_b": label_b,
        "n_clones_a": len(a),
        "n_clones_b": len(b),
        "n_shared_clones": len(shared),
        "jaccard_presence": len(shared) / len(union) if len(union) else np.nan,
        "morisita_horn_abundance": morisita_horn(a, b),
        "fraction_templates_a_in_shared_clones": float(a.loc[shared].sum() / a.sum()) if len(shared) else 0.0,
        "fraction_templates_b_in_shared_clones": float(b.loc[shared].sum() / b.sum()) if len(shared) else 0.0,
        "top20_clone_overlap_fraction": len(top_a & top_b) / 20.0,
    }


def analyze_cardiac_tcr(
    data_dir: pathlib.Path, results_dir: pathlib.Path, tables_dir: pathlib.Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cells = cardiac_cell_clone_table(data_dir, results_dir)
    patient_metrics = patient_clone_metrics(
        cells, "donor", "clone_key", ["condition"], rarefaction_depth=50
    )
    patient_metrics.to_csv(tables_dir / "Table_S33_cardiac_TCR_patient_metrics.tsv", sep="\t", index=False)

    cardiac_tests = group_tests(
        patient_metrics,
        "condition",
        "myocarditis",
        "control",
        [
            "expanded_cell_fraction_ge2",
            "large_clone_cell_fraction_ge5",
            "normalized_clonality",
            "simpson_concentration",
            "top10_fraction",
            "rarefied_clonality_50",
            "rarefied_top10_fraction_50",
        ],
        "cardiac_myocarditis_vs_control_exploratory_only_two_controls",
        minimum_per_group=2,
    )

    sample_meta = parse_soft_samples(data_dir / "GSE228597_family.soft.gz", "GSM8019")
    repertoires: dict[tuple[str, str, str], pd.Series] = {}
    sample_rows: list[dict] = []
    for path in sorted(data_dir.glob("GSM8019*.tsv.gz")):
        gsm = path.name.split("_", 1)[0]
        meta_match = sample_meta.loc[sample_meta["gsm"] == gsm]
        if meta_match.empty:
            continue
        meta = meta_match.iloc[0].to_dict()
        rep = pd.read_csv(path, sep="\t")
        valid = (
            rep["aminoAcid"].notna()
            & (rep["sequenceStatus"].astype(str) == "In")
            & ~rep["aminoAcid"].astype(str).str.contains("\\*", regex=True)
        )
        rep = rep.loc[valid, ["nucleotide", "aminoAcid", "count (templates/reads)"]].copy()
        tissue = str(meta.get("tissue", meta.get("source_name", "unknown")))
        tissue_class = (
            "tumor"
            if tissue.lower() == "tumor"
            else "heart"
            if "heart" in tissue.lower()
            else "adjacent_normal"
            if "normal" in tissue.lower()
            else tissue.lower().replace(" ", "_")
        )
        donor = str(meta.get("donor", "unknown"))
        for clone_definition, clone_column in (
            ("nucleotide_exact", "nucleotide"),
            ("amino_acid_exact_sensitivity", "aminoAcid"),
        ):
            counts = rep.groupby(clone_column)["count (templates/reads)"].sum()
            repertoires[(donor, tissue_class, clone_definition)] = counts
            sample_rows.append(
                {
                    "gsm": gsm,
                    "donor": donor,
                    "tissue": tissue,
                    "tissue_class": tissue_class,
                    "clone_definition": clone_definition,
                    "n_clonotypes": len(counts),
                    "n_templates": int(counts.sum()),
                }
            )

    overlap_rows: list[dict] = []
    donors = sorted({key[0] for key in repertoires})
    for donor in donors:
        for clone_definition in ("nucleotide_exact", "amino_acid_exact_sensitivity"):
            if (donor, "heart", clone_definition) not in repertoires:
                continue
            heart = repertoires[(donor, "heart", clone_definition)]
            for comparator in ("tumor", "adjacent_normal"):
                if (donor, comparator, clone_definition) not in repertoires:
                    continue
                row = {"donor": donor, "clone_definition": clone_definition}
                row.update(
                    overlap_metrics(
                        heart,
                        repertoires[(donor, comparator, clone_definition)],
                        "heart",
                        comparator,
                    )
                )
                overlap_rows.append(row)
    overlap = pd.DataFrame(overlap_rows)
    overlap.to_csv(tables_dir / "Table_S34_heart_tumor_TCR_overlap.tsv", sep="\t", index=False)
    pd.DataFrame(sample_rows).to_csv(
        results_dir / "GSE228597_bulk_TCR_sample_audit.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    paired_rows: list[dict] = []
    if len(overlap):
        pivot_metrics = [
            "jaccard_presence",
            "morisita_horn_abundance",
            "fraction_templates_a_in_shared_clones",
            "top20_clone_overlap_fraction",
        ]
        for clone_definition, definition_block in overlap.groupby("clone_definition", observed=True):
            for metric in pivot_metrics:
                wide = definition_block.pivot(index="donor", columns="compartment_b", values=metric).dropna()
                if not {"tumor", "adjacent_normal"}.issubset(wide.columns) or len(wide) < 3:
                    continue
                try:
                    test = stats.wilcoxon(
                        wide["tumor"], wide["adjacent_normal"], alternative="two-sided", method="auto"
                    )
                    p = float(test.pvalue)
                except ValueError:
                    p = 1.0
                paired_rows.append(
                    {
                        "analysis": f"heart_tumor_vs_heart_adjacent_normal_paired_{clone_definition}",
                        "metric": metric,
                        "group1": "heart_tumor",
                        "group0": "heart_adjacent_normal",
                        "n_group1": len(wide),
                        "n_group0": len(wide),
                        "median_group1": float(wide["tumor"].median()),
                        "median_group0": float(wide["adjacent_normal"].median()),
                        "median_difference": float((wide["tumor"] - wide["adjacent_normal"]).median()),
                        "cliffs_delta": np.nan,
                        "p": p,
                        "clone_definition": clone_definition,
                    }
                )
    paired_tests = pd.DataFrame(paired_rows)
    if len(paired_tests):
        paired_tests["fdr_within_analysis"] = bh_adjust(paired_tests["p"])
    return cardiac_tests, paired_tests


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("project_dir")
    args = parser.parse_args()
    root = pathlib.Path(args.project_dir).resolve()
    data_dir = root / "01_data"
    results_dir = root / "02_results"
    tables_dir = root / "tables_for_article"

    colon_overall, colon_states = analyze_colon_tcr(results_dir, tables_dir)
    coupling_correlations = couple_colon_tcr_to_epithelial_modules(root, results_dir, tables_dir)
    cardiac_overall, cardiac_overlap = analyze_cardiac_tcr(data_dir, results_dir, tables_dir)
    tests = pd.concat(
        [table for table in [colon_overall, colon_states, cardiac_overall, cardiac_overlap] if len(table)],
        ignore_index=True,
        sort=False,
    )
    tests.to_csv(tables_dir / "Table_S35_TCR_statistical_summary.tsv", sep="\t", index=False)
    tests.to_csv(
        results_dir / "TCR_statistical_summary.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    print("TCR evidence analysis complete")
    print(tests.to_string(index=False))
    print("\nWithin-colitis TCR-epithelial coupling (exploratory):")
    print(coupling_correlations.to_string(index=False))


if __name__ == "__main__":
    main()
