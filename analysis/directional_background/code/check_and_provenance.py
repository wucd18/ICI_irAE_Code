from common import *
import subprocess,ast,platform,shutil
from scipy.stats import rankdata
checks=[]
def check(name,ok,detail):
 checks.append(dict(check=name,passed=bool(ok),detail=detail))
 if not ok:raise AssertionError((name,detail))
cl=read(ROOT/'results/classification_stability.tsv');inv=read(BASE/'inputs/comparison_inventory.tsv');planned=inv[inv.condition_class=='ligand']
check('all_planned_keys_once',not cl.duplicated(['ligand','family']).any() and set(zip(cl.ligand,cl.family))==set(zip(planned.ligand,planned.celltype)),len(cl))
check('eligibility_unchanged',np.array_equal(cl.eligible.to_numpy(),planned.eligibility.eq('ELIGIBLE').to_numpy()),int(cl.eligible.sum()))
check('NE_reasons_preserved',cl.loc[~cl.eligible,'exclusion_reason'].notna().all() and (cl.loc[~cl.eligible,['raw_class','bg_relative_class','nonrank_class']]=='NE').all().all(),int((~cl.eligible).sum()))
e=read(ROOT/'results/three_method_donor_effects.tsv');genes,full,sets=definitions()
for r in cl[cl.eligible].itertuples():
 b=e[(e.ligand==r.ligand)&(e.celltype==r.family)]
 for me in ['raw','bg_relative','nonrank']:
  w=b.pivot(index='donor',columns='module',values=me).reindex(DONORS);actual='BOTH' if ((w[UP]<0)&(w[DEP]>0)).sum()==len(DONORS) else 'UP_ONLY' if (w[UP]<0).sum()==len(DONORS) else 'DEPLETED_ONLY' if (w[DEP]>0).sum()==len(DONORS) else 'OTHER_OBSERVED'
  assert actual==getattr(r,me+'_class'),(r.ligand,r.family,me)
check('independent_direction_classification',True,'Recomputed from joint donor rows for each method; strict signs, no rounding')
check('same_donor_pairing_complete',not e.duplicated(['ligand','celltype','donor','module']).any() and e.groupby(['ligand','celltype','module']).donor.nunique().eq(3).all(),len(e))
check('background_subtraction_identity',np.allclose(e.bg_relative,e.raw-e.B,rtol=0,atol=2e-17),float(abs(e.bg_relative-(e.raw-e.B)).max()))
check('new_scores_no_inference_P',not any(c.lower() in ['p','pval','pvalue','q','fdr'] for c in e.columns),'Gene-set reference percentiles retained, no inferential columns')
orig=read(BASE/'results/response_pattern_summary.tsv');o=orig[orig.condition_class=='ligand'].set_index(['ligand','celltype']);check('raw_classes_preserved',all(r.raw_class==o.loc[(r.ligand,r.family),'pattern'] for r in cl.itertuples()),int(cl.original_both.sum()))
check('original_GSEA_byte_preserved',sha(BASE/'results/all_ligand_programme_responses.tsv')==sha(ROOT/'results/original_GSEA_results.tsv'),'No NES/P/q/model refit')
sf=read(ROOT/'review/paired_input_and_normalization_checks.tsv');check('saved_counts_and_factors',sf.pool_counts_equal_saved.all() and sf.saved_factor_extrema_max_abs_error.max()<1e-11,dict(n=len(sf),max_error=sf.saved_factor_extrema_max_abs_error.max()))
libs=read(ROOT/'results/matched_library_index.tsv');done=set()
for r in libs.itertuples():
 if r.status!='ESTIMATED' or r.library_id in done:continue
 z=np.load(ROOT/r.path);a=z['indices'];t=z['target_indices'];assert a.shape==(1000,len(t)) and not np.isin(a,t).any() and all(len(np.unique(v))==len(t) for v in a);done.add(r.library_id)
check('random_libraries_unique_exclude_target',True,dict(unique_libraries=len(done),n_sets_each=1000))
strata=read(ROOT/'results/matching_strata.tsv');check('matching_stratum_capacity',bool((strata.n_candidates>=strata.n_target).all()),'All strata support sampling without replacement; continuous imbalance disclosed separately')
member=read(ROOT/'results/member_gene_log_changes.tsv');check('member_gene_values_formula',np.allclose(member.gene_log_delta,np.log2(member.treatment_normalized_count+1)-np.log2(member.control_normalized_count+1),atol=2e-14,rtol=0),len(member))
ag=member.groupby(['ligand','celltype','donor','module']).gene_log_delta.mean();n=read(ROOT/'results/nonrank_donor_effects.tsv').set_index(['ligand','celltype','donor','module']);check('nonrank_unweighted_all_members',np.allclose(ag,n.loc[ag.index,'normalized_expression_change'],atol=2e-14,rtol=0),'No significance-gene filtering')
counts=member.groupby(['ligand','celltype','donor','module']).gene.nunique();check('fixed_member_coverage',all(int(count)==len(sets[key[-1]]) for key,count in counts.items()),{m:len(ix) for m,ix in sets.items()});del member
sub=read(ROOT/'results/within_family_substate_proportions.tsv');supported=sub[sub.common_supported];w=supported.groupby(['ligand','celltype','donor']).reference_weight.sum();check('composition_fixed_weights_sum1',np.allclose(w,1,rtol=0,atol=1e-14),len(w))
ce=read(ROOT/'results/composition_standardized_expression.tsv');valid=ce[ce.status=='ESTIMATED'];check('composition_support_requirement',bool((valid.treatment_retained_fraction>=.8).all() and (valid.control_retained_fraction>=.8).all()),len(valid));check('composition_NE_not_imputed',ce.loc[ce.status!='ESTIMATED','standardized_common_delta'].isna().all(),int((ce.status!='ESTIMATED').sum()))
pr=read(ROOT/'results/pool_relative_effects.tsv');check('pool_reference_separate',pr.loc[pr.status=='ESTIMATED','n_reference_ligands'].ge(3).all() and not pr.reference_ligands.str.contains('No_Cytomix').any(),'No background+pool double subtraction')
used=[TASK]
for p in [BASE/'inputs/frozen_module_definitions.csv',BASE/'inputs/IBD_membership_definitions.tsv',BASE/'inputs/comparison_inventory.tsv',BASE/'inputs/all_ligand_pool_eligibility.tsv',BASE/'inputs/all_ligand_pseudobulk_metadata.tsv',BASE/'inputs/pseudobulk_genes.tsv',BASE/'inputs/new_model_jobs.tsv',BASE/'inputs/targeted_model_inventory.tsv',BASE/'results/new_model_inventory.tsv',BASE/'results/all_readout_donor_deltas_long.tsv',BASE/'results/response_pattern_summary.tsv',BASE/'results/all_ligand_programme_responses.tsv',BASE/'code/stream_r_object.py',BASE/'code/prepare_increment.py',BASE/'code/score_ibd_references.py',BASE/'code/increment_models.R',EXT/'source_data/pseudobulk_counts.mtx',EXT/'source_data/pseudobulk_metadata.tsv',RAW/'counts_parsed.pkl',RAW/'cell_metadata.tsv.gz',RAW/'column_pointers.npy',RAW/'genes.tsv',OLD/'provenance/Capeling_Figure_4_Code_pinned.Rmd',OLD/'provenance/Cytomix_Secretome_Processing.Rmd',OLD/'provenance/Capeling_publisher.html',Path(r'E:\ICI_irAE_Original_Project\01_data\GSE313368_article_supplementary.zip')]:used.append(p)
jobs=read(BASE/'inputs/new_model_jobs.tsv')
for r in jobs.itertuples():used.extend([Path(r.counts),Path(r.metadata)])
approved=read(BASE/'review/input_manifest.tsv');oldmap=approved.set_index('path').to_dict('index');manifest=read(BASE/'SHA256_manifest.tsv');basemap=manifest.set_index('relative_path').to_dict('index')
for p in approved.path:
 if p.endswith('_counts_i.bin') or p.endswith('_counts_x.bin'):used.append(Path(p))
records=[]
for p in sorted(set(used),key=str):
 assert p.is_file(),p;size=p.stat().st_size;prior=oldmap.get(str(p));relative=p.relative_to(BASE).as_posix() if p.is_relative_to(BASE) else None;prior=basemap.get(relative,prior)
 reusebig=size>1_000_000_000 and prior is not None
 h=prior['sha256'] if reusebig else sha(p)
 if prior is not None:assert size==prior['bytes'] and (reusebig or h==prior['sha256']),(p,'changed approved input')
 records.append(dict(path=str(p),bytes=size,sha256=h,source='GSE313368 / approved all-ligand run' if p!=TASK else 'User task file',role='matching/counts/annotation/reference/scoring or approved code' if p!=TASK else 'authorized scope',hash_verification='REUSED_APPROVED_SHA_SIZE_CHECK_ONLY_NOT_REHASHED_THIS_ROUND' if reusebig else 'CURRENT_SHA256_COMPUTED',approved_manifest_match='SIZE_ONLY_INHERITED_SHA' if reusebig else 'MATCH' if prior is not None else 'NEWLY_RECORDED'))
save(pd.DataFrame(records),'review/input_manifest.tsv');check('scoped_input_identity',True,dict(inputs=len(records),large_blob_hashes_inherited=sum(r['hash_verification'].startswith('REUSED') for r in records)))
# Syntax and prohibited operations are checked only in new analysis code, not the copied legacy parser.
static=[]
for p in (ROOT/'code').glob('*.py'):
 ast.parse(p.read_text(encoding='utf8'));static.append(dict(file=p.name,syntax='PASS',sha256=sha(p)))
save(pd.DataFrame(static),'review/static_source_checks.tsv');check('python_static_syntax',True,len(static))
check('text_outputs_present',all((ROOT/'text'/p).is_file() for p in ['results_update.md','figure_legends_update.md','conclusions_update.md','claim_evidence_matrix.tsv']),'Original manuscript unmodified; standalone proposals')
save(pd.DataFrame(checks),'review/targeted_numerical_checks.tsv');jw(ROOT/'review/targeted_validation.json',dict(status='PASS',scope='Only new directional/background sensitivity products; not full-project audit',checks=checks,completed=now(),environment=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__)))
jw(ROOT/'review/KNOWN_EXECUTION_NOTES.json',dict(initial_fixture='One initial fixture expected floating-point standard deviation to equal mathematical zero with atol=0; observed <3e-14. Fixture tolerance corrected to 1e-12. Scientific sign tests were not changed. Original failed log preserved.',empty_merge_table='No sparse joint bins required merging. Empty table now has its declared header; no random library or result changed.',analysis_model_refits=False,new_inferential_tests=False,source_hash_scope='Large original count blobs inherit approved SHA-256 with current size checks. All smaller used inputs and approved baseline counts/results were currently hashed; no claim of a new full source audit.'))
print('TARGETED_VALIDATION_PASS',len(checks),'checks',len(records),'inputs',flush=True)
