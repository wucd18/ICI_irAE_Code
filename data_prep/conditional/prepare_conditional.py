import os
from pathlib import Path
import importlib.util,hashlib,json,sys
import numpy as np,pandas as pd
from datetime import datetime,timezone
O=Path(__file__).resolve().parents[1];P=Path(os.environ['ICI_ARCHIVE_ROOT'])/'original_v1'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(4*1024*1024):h.update(b)
 return h.hexdigest()
paths=[Path(os.environ['ICI_ORIGINAL_CODE_ROOT'])/'00_scripts/evaluate_patient_level_module_performance.py',P/'tables_for_article/Table_S21_frozen_module_definitions.tsv',P/'tables_for_article/Table_S13_GSE206300_cluster_family_mapping.tsv',P/'tables_for_article/Table_S55_patient_level_frozen_module_scores.tsv']
for ds in ['GSE206300','GSE253720']:
 for suffix in ['metadata.tsv','genes.tsv','counts.mtx.gz']:paths.append(P/'02_results'/f'{ds}_pseudobulk_{suffix}')
manifest=[dict(path=str(p),size=p.stat().st_size,sha256=sha(p),role='immutable V1 input for conditional extension') for p in paths]
plan=dict(utc=datetime.now(timezone.utc).isoformat(),authorization='User explicitly requested reanalysis and changed strategy in current turn.',primary_dataset='GSE206300',comparison='PD-1 only irColitis vs On ICI therapy',unit='patient within the four original epithelial families; not four independent cohorts',outcomes=['ICI_COLITIS_EPITHELIAL_UP_NON_IFN','ICI_COLITIS_EPITHELIAL_LOSS_STRICT'],scores='Original frozen definitions, mean within-patient gene rank; loss oriented as 1 minus raw rank.',model='oriented score ~ disease indicator + centered IFN score',inference='OLS coefficients; HC3 standard errors and two-sided residual-df t approximation; BH over all 8 primary family-by-module tests. 95% intervals are pointwise.',discovery='GSE253720 same two outcomes; separately adjusted self-check only, excluded from independent validation.',sensitivity='Leave each patient out once, retain every coefficient, report sign stability without additional significance claims.',no_changes='No gene re-selection, candidate re-ranking, threshold changes, or clinical diagnostic/safety claims.',limitations='Small sample observational conditional associations. IFN may be a disease mediator; adjustment is descriptive, not causal independence. Collinearity and IFN score overlap must be reported.',inputs=manifest)
lp=O/'review/conditional_analysis_plan_locked.json'
with lp.open('x',encoding='utf-8') as f:json.dump(plan,f,indent=2)
pd.DataFrame(manifest).to_csv(O/'review/conditional_inputs_before.tsv',sep='\t',index=False)
spec=importlib.util.spec_from_file_location('v1_score_functions',paths[0]);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
defs=pd.read_csv(paths[1],sep='\t');mapping=pd.read_csv(paths[2],sep='\t')
all_s=[]
for ds in ['GSE206300','GSE253720']:
 m,meta,genes=mod.read_pb(P,ds+'_pseudobulk')
 if ds=='GSE206300':m,meta=mod.aggregate_independent_families(m,meta,mapping);role='independent_primary_PD1_only'
 else:
  keep=meta.celltype.eq('Epithelial')&meta.condition.isin(['CPI_colitis','HC']);m=m[:,np.flatnonzero(keep)];meta=meta.loc[keep].copy();meta['sample_id']=meta.pseudobulk_id;role='module_derivation_self_check'
 assert not meta.duplicated(['patient_id','celltype']).any()
 assert meta.n_cells.ge(20).all()
 scores=mod.rank_scores(m,genes,meta,defs,ds,role);all_s.append(scores)
s=pd.concat(all_s,ignore_index=True)
original=pd.read_csv(paths[3],sep='\t');keys=['dataset','sample_id','module']
c=s.merge(original,on=keys,suffixes=('_new','_old'),validate='one_to_one')
assert len(c)==len(s)
for col in ['raw_mean_within_sample_gene_rank','disease_oriented_score','n_genes_available']:
 assert np.allclose(c[col+'_new'],c[col+'_old'],rtol=0,atol=1e-14),col
s.to_csv(O/'source_data/conditional_patient_scores.tsv',sep='\t',index=False)
wide=s.pivot(index=['dataset','cohort_role','patient_id','sample_id','condition','celltype','n_cells_or_spots'],columns='module',values='disease_oriented_score').reset_index()
assert wide.notna().all().all()
wide.to_csv(O/'source_data/conditional_model_input.tsv',sep='\t',index=False)
(O/'review/conditional_score_reproduction.json').write_text(json.dumps(dict(status='PASS',rows=len(s),tolerance=1e-14,columns=3,all_input_hashes_unchanged=all(sha(Path(x['path']))==x['sha256'] for x in manifest)),indent=2),encoding='utf-8')
print('FROZEN_SCORE_REPRODUCTION_PASS',len(s))
