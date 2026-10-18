from common import *
import zipfile,platform,shutil
cfg=jr(ROOT/'analysis_config.json');g=pd.read_csv(ROOT/'results/all_ligand_programme_responses.tsv',sep='\t');summary=pd.read_csv(ROOT/'results/response_pattern_summary.tsv',sep='\t');new=pd.read_csv(ROOT/'review/new_model_checks_and_provenance.tsv',sep='\t')
assert jr(ROOT/'review/product_numerical_checks.json')['status']=='PASS'
assert jr(ROOT/'review/visual_review.json')['status']=='PASS'
logs=[jr(p) for p in sorted((ROOT/'logs').glob('execution_enrichment_*.json'))];assert len(logs)==4 and all(x['completed_records_exist'] for x in logs)
runtime_ok=all(x['exit_code']==0 for x in logs)
status=dict(status='COMPLETE_WITH_DISCLOSED_LIMITATIONS',analysis_outputs='PASS_NEW_OUTPUT_CHECKS',visual_outputs='PASS_RENDERED_AND_INSPECTED',R_process_exit='PASS' if runtime_ok else 'NONZERO_EXIT_AFTER_COMPLETE_OUTPUTS_NOT_PASS',new_gene_models=len(new),reused_targeted_models=25,reused_targeted_module_results=100,primary_planned=948,primary_estimable=834,primary_NE=114,ligand_comparisons_planned=316,ligand_comparisons_eligible=278,ligand_comparisons_NE=38,pdf_figures=6,pdf_pages=15,manuscript_original_unchanged=True,full_manuscript_overwritten=False,public_upload=False,analysis_incomplete_items=[],unresolved_items=[dict(id='IL22_DOSE_METADATA',status='BLOCKED_SOURCE_CONFLICT',impact='No IL22 concentration assigned; published dose table lists IL11 whereas released data/code label IL22. Analysis follows the released label and remains in atlas.'),dict(id='R_RUNTIME_EXIT',status='UNRESOLVED_ENVIRONMENT_ISSUE',impact='R workers emitted completion records and complete validated outputs, then returned 3221225477 (0xC0000005). Base R exits 0; loading fgsea alone reproduces the exit fault. Underlying cause not established. No package installation or refitting performed.')],recovery=dict(failed_run=str(ROOT.with_name(ROOT.name.replace('_r3','_r2'))),cause='CSV BOM header parsing yielded empty gene sets',resolution='New run explicitly decodes UTF-8-BOM, checks every member identity, reuses unchanged new gene statistics, recomputes only enrichment. Invalid r2 enrichment excluded.'),created_utc=now())
jw(ROOT/'review/FINAL_STATUS.json',status)
jw(ROOT/'review/python_environment.json',dict(python=sys.version,executable=sys.executable,platform=platform.platform(),numpy=np.__version__,pandas=pd.__version__))
archive=PROJECT/(ROOT.name+'_DELIVERY.zip')
readme=f'''# CSBJ all-ligand incremental analysis delivery

## Main result and status

The requested incremental analysis, six figures, source tables, Results/Methods/figure-legend/Conclusions patches and Chinese analysis report are complete. The package is **COMPLETE_WITH_DISCLOSED_LIMITATIONS**, not an unqualified software-runtime PASS or automatic journal submission.

Run: `{ROOT}`

Archive: `{archive}`

Start with `text/ANALYSIS_REPORT_ZH.md`, then `text/Results_additions.md` and `figures/Figure4_TNFSF12.pdf`.

## What was actually computed

- All 79 released ligands were retained. There are 316 planned ligand–family combinations, 278 eligible and 38 NE; 78 ligands have at least one eligible family. No-cytomix is separate.
- Twenty-one approved ligand models plus four assay-reference models were reused. Only 257 missing eligible gene models were fitted, once each. No patient, TCR, tumour, AUC or IFN-conditional analysis was rerun.
- Primary BH: 948 planned tests, 834 estimated, 114 NE. Residual BH: 316 planned, 278 estimated, 38 NE. Assay-reference BH: 16. Old targeted FDR stays in separate columns.
- Unrounded same-donor directions: 8 BOTH, 87 UP_ONLY, 7 DEPLETED_ONLY, 176 OTHER_OBSERVED among 278 eligible comparisons. Zero never counts as positive or negative. OTHER does not mean no response.
- Same-rule IFN and published IBD membership benchmarks include mean and donor-stratified descriptive correlations without independence-based P values. IBD membership rank scores are explicitly not the publisher's original activity/cNMF scores.

## Required deliverables

| File | Contents |
| --- | --- |
| results/all_ligand_programme_responses.tsv | Every planned ligand/reference × family × programme; NES, P, atlas FDR, old targeted FDR, score summaries, coverage, dose, origin |
| results/all_ligand_donor_deltas.tsv | Planned donor rows, actual paired scores/directions and missingness; raw depleted delta is not inverted |
| results/response_pattern_summary.tsv | Same-donor joint categories, sign counts, means, medians, ranges and denominators |
| results/programme_benchmark.tsv | Within-family Pearson/Spearman and common-n, including each donor stratum |
| results/non_estimable_contrasts.tsv | All 38 NE ligand–family combinations and donor-coverage reason |
| figures/ | Six final PDFs and editable per-page SVGs |
| figure_source_data/ | Exact plot inputs and prespecified example-selection assignments |
| text/ | Results, Methods, legends, Conclusions options, change map and Chinese analysis report |
| review/ | Input hashes, model checks, numerical and visual checks, runtime diagnosis, recovery record |
| logs/ | Original model-fit and enrichment-recovery commands, times, stdout, stderr, exit codes and R sessionInfo |

## Interpretation and unresolved items

1. **IL22 dose:** published labels/code say IL22; the source dose table instead lists IL11. Its dose is unresolved, not inferred. The 78 other protein matches and every original dose-table entry are preserved. Vendor/catalogue text is copied verbatim; apparent source typographical inconsistencies were not silently corrected.
2. **R process exit:** each analysis worker wrote all completion records and all output tables, then returned Windows access-violation code 0xC0000005. New output/schema checks, independent BH verification and input-immutability checks pass. Process exits do not pass. Minimal diagnostics reproduce the fault when fgsea is merely loaded; its root cause remains unresolved. This should be fixed in a separately authorized environment-maintenance step before claiming a clean-runtime release. No package was installed or upgraded here.
3. **IBD definition:** CA1 occurs twice in the published IBD2 member list and is counted once. Exact-symbol coverage is 106/113 and 86/91. Missing symbols are listed, not silently renamed. Weighted/activity scores were not pretended to have been reproduced.
4. **Scope:** these are responses under the tested donor/pool/cell-family conditions, not potency, protection, clinical efficacy, family-interaction or safety evidence. Existing negative, discordant, NE and tumour/off-target guardrails remain applicable.

## Recovery without refitting

The first preparation attempt stopped at a pickle import before fitting. Run r2 then completed the 257 gene models, but its R CSV reader did not decode the BOM and therefore passed empty gene sets to enrichment. That was a technical failure, not biological NE. Run r3 preserves/reuses those exact gene-statistic bytes, explicitly decodes UTF-8-BOM, checks all member identities and regenerates only the enrichment stage. All invalid r2 GSEA outputs are excluded from final result tables. The failed run is retained separately; representative error evidence and original logs are included.

## Current manuscript and historical Figure 4c

`inputs/Manuscript_CSBJ_revised_clean.docx` is an unchanged copy of the current CSBJ manuscript. No old Word/Markdown overwrites it. Its original Figure 4 image, legend and selected-set context are preserved in `inputs/historical_Figure4/`. The delivered text files contain actual ready-to-use additions/replacements based on results; they are not a newly assembled full Word manuscript. Descriptively named new figures await the author's main/supplementary placement choice; unchanged figures have not been silently renumbered. The old selected-set transfer is proposed for supplementary placement, with original selection caveats intact.

## Reuse and reproducibility

Use the existing result tables for plotting and text (`plot_increment.py`, `write_manuscript_additions.py`). `summarize_scores.py` and `combine_gsea.py` derive the final response/benchmark tables from saved results. `recover_enrichment.R` documents exactly the executed enrichment repair and refuses to overwrite completed results. `increment_models.R` contains the repaired member-parser checks for a future newly isolated run; it was not used to refit models here. Do not run RUN_ALL or restart the full pipeline to inspect these results.

The archive includes all newly aggregated model count matrices, all 257 saved gene-statistic tables, fixed member definitions, reused targeted summaries/scores, code, figures and QA evidence. It does not duplicate the multi-gigabyte original serialized/sparse source dataset; original paths and SHA-256 hashes are in `review/input_manifest.tsv` and the prior locked plan. Regenerating missing primary source data is outside this incremental delivery. Existing approved scientific executables were used; bundled Python/Poppler was used only for publisher workbook reading and PDF QA.

## Integrity

`SHA256_manifest.tsv` lists every delivered file except itself. Temporary R folders and Python bytecode are excluded. The external ZIP verification receipt records CRC verification, every archived member hash against this manifest, and the archive SHA-256. No upload, GitHub push, Figshare release or journal submission was performed.
'''
(ROOT/'README_DELIVERY.md').write_text(readme,encoding='utf8')
files=[p for p in ROOT.rglob('*') if p.is_file() and not any(part in ['temp','__pycache__'] for part in p.relative_to(ROOT).parts) and p.name!='SHA256_manifest.tsv']
records=[dict(relative_path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(files)]
save(pd.DataFrame(records),'SHA256_manifest.tsv')
files.append(ROOT/'SHA256_manifest.tsv')
assert not archive.exists(),'Refusing to overwrite an existing delivery archive'
with zipfile.ZipFile(archive,'w',allowZip64=True) as z:
    for i,p in enumerate(sorted(files)):
        typ=zipfile.ZIP_STORED if p.suffix.lower() in ['.gz','.png','.docx','.xlsx'] else zipfile.ZIP_DEFLATED
        z.write(p,arcname=ROOT.name+'/'+p.relative_to(ROOT).as_posix(),compress_type=typ,compresslevel=4 if typ==zipfile.ZIP_DEFLATED else None)
        if i%300==0:print('PACKAGED',i,'/',len(files),flush=True)
manifest={r['relative_path']:r for r in records}
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    expected={ROOT.name+'/'+n for n in manifest}|{ROOT.name+'/SHA256_manifest.tsv'};assert set(z.namelist())==expected
    for n,r in manifest.items():assert hashlib.sha256(z.read(ROOT.name+'/'+n)).hexdigest()==r['sha256'],n
receipt=dict(status='PASS_ARCHIVE_INTEGRITY',archive=str(archive),bytes=archive.stat().st_size,sha256=sha(archive),file_count=len(files),all_manifest_members_match=True,zip_crc_pass=True,scientific_and_runtime_status='See included review/FINAL_STATUS.json; integrity PASS does not mean R process-exit PASS',created_utc=now())
jw(archive.with_suffix('.verification.json'),receipt)
archive.with_suffix('.sha256').write_text(receipt['sha256']+'  '+archive.name+'\n',encoding='ascii')
print('DELIVERY_COMPLETE',json.dumps(receipt),flush=True)
