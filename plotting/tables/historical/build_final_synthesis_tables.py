#!/usr/bin/env python3
"""Build manuscript-facing cohort, gate, candidate, and claim-evidence tables."""

from pathlib import Path
import argparse
import json
import os
import numpy as np
import pandas as pd


def truth(value):
    if pd.isna(value):
        return False
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    return bool(value)


def yn(value):
    return "yes" if truth(value) else "no"


def load_annotations(path):
    def unique_object(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError('Duplicate annotation key: ' + key)
            result[key] = value
        return result
    annotations = json.loads(Path(path).read_text(encoding='utf-8-sig'), object_pairs_hook=unique_object)
    required = {'ligands', 'gate_rows', 'claims', 'candidate_decisions'}
    if not isinstance(annotations, dict) or not required.issubset(annotations):
        raise ValueError('BLOCKED_INPUT: incomplete versioned historical annotations')
    ligands = annotations['ligands']
    if not isinstance(ligands, list) or not ligands or not all(isinstance(x, str) and x for x in ligands) or len(set(ligands)) != len(ligands):
        raise ValueError('Missing/duplicate archived candidate list')
    for key in ('gate_rows', 'claims'):
        if not isinstance(annotations[key], list) or not annotations[key] or any(not isinstance(x, list) or len(x) != 5 for x in annotations[key]):
            raise ValueError('Malformed archived annotation rows: ' + key)
    decisions = annotations['candidate_decisions']
    if not isinstance(decisions, dict) or not all(isinstance(decisions.get(k), str) and decisions[k] for k in ligands + ['SORL1']):
        raise ValueError('Missing archived candidate decision')
    return annotations


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("project_dir", type=Path)
    args = parser.parse_args()
    root = args.project_dir.resolve()
    tables = root / "tables_for_article"

    annotation_path = Path(os.environ['ICI_HISTORICAL_ANNOTATIONS'])
    annotations = load_annotations(annotation_path)
    catalog = pd.read_csv(tables / "Table_S1_dataset_catalog.tsv", sep="\t")
    used = catalog[catalog["final_role"] != "not used"].copy()
    table1 = used[[
        "accession", "organ_or_context", "condition", "assays", "analysis_n",
        "control_design", "final_role", "decision_note",
    ]].rename(columns={
        "accession": "cohort",
        "organ_or_context": "context",
        "analysis_n": "evaluable_samples",
        "final_role": "role_in_study",
        "decision_note": "key_constraint",
    })
    table1.to_csv(tables / "Table_1_study_cohorts_and_roles.tsv", sep="\t", index=False)

    screen = pd.read_csv(tables / "Table_S23_GSE313368_irAE_rescue_screen.tsv", sep="\t")
    guard = pd.read_csv(tables / "Table_S30_candidate_tumor_efficacy_guardrail_summary.tsv", sep="\t")
    off = pd.read_csv(tables / "Table_S48_organoid_candidate_offtarget_summary.tsv", sep="\t")
    downstream = pd.read_csv(tables / "Table_S43_TNFSF12_downstream_multicohort_ranking.tsv", sep="\t")

    ligands = annotations["ligands"]
    rows = []
    for candidate in ligands:
        s = screen.loc[screen["perturbation"] == candidate].iloc[0]
        g = guard.loc[guard["candidate"] == candidate].iloc[0]
        o = off.loc[off["perturbation"] == candidate].iloc[0]
        inflammatory = int(o["n_positive_inflammatory_fdr05"]) > 0
        ifn_loss = int(o["n_negative_IFN_fdr05"]) > 0
        tumor_fail = str(g["guardrail_status"]).startswith("FAIL")
        expanded_pass = (
            truth(s["function_rescue_gate"])
            and truth(s["disease_program_not_aggravated_gate"])
            and not inflammatory
            and not ifn_loss
            and str(g["guardrail_status"]) == "NO_DETECTED_HARM_not_equivalent_to_safety"
        )

        reasons = []
        if not truth(s["disease_program_not_aggravated_gate"]):
            reasons.append("aggravates a prespecified disease program")
        if inflammatory:
            reasons.append("activates inflammatory Hallmark pathways")
        if ifn_loss:
            reasons.append("suppresses at least one interferon Hallmark pathway")
        if tumor_fail:
            reasons.append("detected adverse tumour-ICI association")
        elif str(g["guardrail_status"]).startswith("CONCERN"):
            reasons.append("tumour guardrail concern")
        if not reasons:
            reasons.append("absence of detected efficacy harm is not evidence of safety")

        rows.append({
            "candidate": candidate,
            "candidate_level": "systemic ligand perturbation",
            "organoid_function_rescue": yn(s["function_rescue_gate"]),
            "legacy_module_decoupling_gate": yn(s["full_decoupling_gate"]),
            "expanded_offtarget_gate": "pass" if expanded_pass else "fail",
            "positive_inflammatory_pathways_fdr05": int(o["n_positive_inflammatory_fdr05"]),
            "negative_interferon_pathways_fdr05": int(o["n_negative_IFN_fdr05"]),
            "tumor_efficacy_guardrail": g["guardrail_status"],
            "independent_tissue_support": "not applicable to the administered ligand as a whole",
            "final_decision": annotations["candidate_decisions"][candidate],
            "decision_basis": "; ".join(reasons),
            "claim_allowed": "mechanistic clue/local downstream program only",
        })

    sorl1 = downstream.loc[downstream["gene"] == "SORL1"].iloc[0]
    rows.append({
        "candidate": "SORL1",
        "candidate_level": "TNFSF12 downstream repair-core gene",
        "organoid_function_rescue": (
            f"up log2FC={sorl1['organoid_log2FC']:.3f}, FDR={sorl1['organoid_fdr']:.3g}"
        ),
        "legacy_module_decoupling_gate": "not applicable",
        "expanded_offtarget_gate": "not established for direct manipulation",
        "positive_inflammatory_pathways_fdr05": np.nan,
        "negative_interferon_pathways_fdr05": np.nan,
        "tumor_efficacy_guardrail": (
            f"response g={sorl1['response_meta_g']:.3f} "
            f"[{sorl1['response_ci_low']:.3f},{sorl1['response_ci_high']:.3f}]; "
            f"OS HR={sorl1['os_meta_hr']:.3f} "
            f"[{sorl1['os_ci_low']:.3f},{sorl1['os_ci_high']:.3f}]"
        ),
        "independent_tissue_support": (
            f"down in {int(sorl1['independent_n_down'])}/"
            f"{int(sorl1['independent_n_families'])} epithelial families; "
            f"spatial logFC={sorl1['spatial_logFC']:.3f}, FDR={sorl1['spatial_fdr']:.3g}"
        ),
        "final_decision": annotations["candidate_decisions"]["SORL1"],
        "decision_basis": (
            "convergent organoid and independent-colitis evidence, but no direct SORL1 "
            "perturbation and no established systemic or tumour safety"
        ),
        "claim_allowed": "candidate marker/mechanistic hypothesis; not a therapeutic target",
    })
    final = pd.DataFrame(rows)
    final.to_csv(tables / "Table_2_final_candidate_decisions.tsv", sep="\t", index=False)
    final.to_csv(tables / "Table_S59_final_candidate_decision_matrix.tsv", sep="\t", index=False)

    gate_rows = annotations["gate_rows"]
    gates = pd.DataFrame(gate_rows, columns=["gate_id", "gate", "final_status", "evidence_based_interpretation", "supporting_outputs"])
    gates.to_csv(tables / "Table_S60_final_gate_status.tsv", sep="\t", index=False)

    claims = pd.DataFrame(annotations["claims"], columns=["claim_id", "candidate_manuscript_claim", "status", "supporting_outputs", "mandatory_boundary"])
    claims.to_csv(tables / "Table_S61_claim_evidence_map.tsv", sep="\t", index=False)

    print(f"Wrote {len(table1)} cohort rows, {len(final)} candidate rows, {len(gates)} gates, and {len(claims)} claim rows")


if __name__ == "__main__":
    main()
