#!/usr/bin/env python
"""Cross-organ convergence and divergence in non-immune tissue receivers."""

from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats


RECEIVERS = {
    "GSE228597": ["Endothelial cells", "Fibroblasts", "Mural cells"],
    "GSE253720": ["Endothelial", "Epithelial", "Mesenchymal Stromal"],
}
ORGAN = {"GSE228597": "heart", "GSE253720": "colon"}


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: analyze_cross_organ_receiver_programs.py <project_dir>")
    root = pathlib.Path(sys.argv[1]).resolve()
    tables = root / "tables_for_article"
    enrich = pd.read_csv(
        root / "02_results" / "celltype_pathway_enrichment_all.tsv.gz", sep="\t"
    )
    blocks = []
    for dataset, celltypes in RECEIVERS.items():
        block = enrich.loc[
            enrich["dataset"].eq(dataset)
            & enrich["celltype"].isin(celltypes)
            & enrich["collection"].eq("HALLMARK")
        ].copy()
        block["receiver_organ"] = ORGAN[dataset]
        block["fdr05_direction"] = np.select(
            [block["padj"].le(0.05) & block["NES"].gt(0), block["padj"].le(0.05) & block["NES"].lt(0)],
            ["up", "down"],
            default="not_significant",
        )
        blocks.append(block)
    receiver = pd.concat(blocks, ignore_index=True)
    receiver.to_csv(
        tables / "Table_S52_cross_organ_receiver_Hallmark_matrix.tsv", sep="\t", index=False
    )

    rows = []
    for pathway, block in receiver.groupby("pathway", observed=True):
        row = {"pathway": pathway}
        for organ in ["heart", "colon"]:
            ob = block.loc[block["receiver_organ"].eq(organ)]
            row[f"{organ}_n_receiver_compartments"] = len(ob)
            row[f"{organ}_n_fdr_up"] = int((ob["fdr05_direction"] == "up").sum())
            row[f"{organ}_n_fdr_down"] = int((ob["fdr05_direction"] == "down").sum())
            row[f"{organ}_median_NES"] = float(ob["NES"].median())
            row[f"{organ}_receiver_results"] = ";".join(
                f"{r.celltype}:{r.NES:.3f}:{r.padj:.3g}" for r in ob.itertuples()
            )
        hsign = np.sign(row["heart_median_NES"])
        csign = np.sign(row["colon_median_NES"])
        row["shared_up_in_at_least_two_receivers_per_organ"] = (
            row["heart_n_fdr_up"] >= 2 and row["colon_n_fdr_up"] >= 2
        )
        row["shared_down_in_at_least_two_receivers_per_organ"] = (
            row["heart_n_fdr_down"] >= 2 and row["colon_n_fdr_down"] >= 2
        )
        row["opposite_organ_median_direction"] = bool(hsign != 0 and csign != 0 and hsign != csign)
        row["cross_organ_class"] = (
            "shared_up_receiver_program"
            if row["shared_up_in_at_least_two_receivers_per_organ"]
            else "shared_down_receiver_program"
            if row["shared_down_in_at_least_two_receivers_per_organ"]
            else "organ_divergent_receiver_program"
            if row["opposite_organ_median_direction"]
            and (row["heart_n_fdr_up"] + row["heart_n_fdr_down"] >= 1)
            and (row["colon_n_fdr_up"] + row["colon_n_fdr_down"] >= 1)
            else "not_classified"
        )
        rows.append(row)
    summary = pd.DataFrame(rows)
    class_order = {
        "shared_up_receiver_program": 1,
        "shared_down_receiver_program": 2,
        "organ_divergent_receiver_program": 3,
        "not_classified": 4,
    }
    summary["class_order"] = summary["cross_organ_class"].map(class_order)
    summary = summary.sort_values(
        ["class_order", "heart_n_fdr_up", "colon_n_fdr_up", "pathway"],
        ascending=[True, False, False, True],
    )
    summary.to_csv(
        tables / "Table_S53_cross_organ_receiver_program_summary.tsv", sep="\t", index=False
    )

    valid = summary[["heart_median_NES", "colon_median_NES"]].dropna()
    corr = stats.spearmanr(valid["heart_median_NES"], valid["colon_median_NES"])
    audit = pd.DataFrame(
        [
            {
                "analysis": "cross_organ_receiver_Hallmark_median_NES",
                "n_pathways": len(valid),
                "spearman_rho": float(corr.statistic),
                "p": float(corr.pvalue),
                "n_shared_up": int((summary["cross_organ_class"] == "shared_up_receiver_program").sum()),
                "n_shared_down": int((summary["cross_organ_class"] == "shared_down_receiver_program").sum()),
                "n_organ_divergent": int(
                    (summary["cross_organ_class"] == "organ_divergent_receiver_program").sum()
                ),
                "interpretation": "shared attack programs coexist with organ-specific receiver-state remodeling",
            }
        ]
    )
    audit.to_csv(
        tables / "Table_S54_cross_organ_receiver_statistical_summary.tsv", sep="\t", index=False
    )
    print("Cross-organ receiver analysis complete")
    print(
        summary.loc[summary["cross_organ_class"].ne("not_classified"),
                    ["pathway", "heart_median_NES", "colon_median_NES", "cross_organ_class"]]
        .to_string(index=False)
    )
    print(audit.to_string(index=False))


if __name__ == "__main__":
    main()
