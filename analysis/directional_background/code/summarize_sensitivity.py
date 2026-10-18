from common import *
import itertools
from scipy.stats import spearmanr
genes,full,sets=definitions();inv=read(BASE/'inputs/comparison_inventory.tsv')
bg=read(ROOT/'results/matched_background_effects.tsv');non=read(ROOT/'results/nonrank_donor_effects.tsv')
keys=['ligand','celltype','donor','module']
effect=bg.rename(columns={'delta_raw':'raw','matched_background_relative_delta':'bg_relative'}).merge(non[keys+['normalized_expression_change']].rename(columns={'normalized_expression_change':'nonrank'}),on=keys,validate='one_to_one')
save(effect,'results/three_method_donor_effects.tsv')
raw=read(ROOT/'results/raw_donor_effects.tsv');raw['control_ids']=raw.apply(lambda r:';'.join('|'.join([r.donor,p,r.celltype,'Control']) for p in str(r.pools).split(';')),axis=1)
raw=raw.merge(inv[['ligand','celltype','eligibility','condition_class']],on=['ligand','celltype'],validate='many_to_one');save(raw,'results/raw_donor_effects.tsv')
def signstats(v):
 v=np.asarray(v,dtype=float);v=v[np.isfinite(v)]
 return dict(n=len(v),mean=float(v.mean()) if len(v) else np.nan,median=float(np.median(v)) if len(v) else np.nan,min=float(v.min()) if len(v) else np.nan,max=float(v.max()) if len(v) else np.nan,n_negative=int((v<0).sum()),n_positive=int((v>0).sum()),n_zero=int((v==0).sum()),negative_fraction=float((v<0).mean()) if len(v) else np.nan,positive_fraction=float((v>0).mean()) if len(v) else np.nan,zero_fraction=float((v==0).mean()) if len(v) else np.nan)
rs=[];ns=[];bs=[]
eligible=effect[effect.condition_class=='ligand']
for groupcols in [['module'],['module','celltype'],['module','celltype','donor'],['module','celltype','donor','pools','control_ids']]:
 for key,b in eligible.groupby(groupcols):
  key=key if isinstance(key,tuple) else (key,);base=dict(zip(groupcols,key));base['stratum']='|'.join(groupcols)
  rs.append(dict(**base,**signstats(b.raw)))
  ns.append(dict(**base,**signstats(b.nonrank)))
  bs.append(dict(**base,**signstats(b.B),random_set_negative_fraction_mean=b.random_negative_fraction.mean(),random_set_positive_fraction_mean=b.random_positive_fraction.mean(),background_relative_negative_fraction=(b.bg_relative<0).mean(),background_relative_positive_fraction=(b.bg_relative>0).mean(),raw_negative_fraction=(b.raw<0).mean(),raw_positive_fraction=(b.raw>0).mean()))
save(pd.DataFrame(rs),'results/raw_direction_summary.tsv');save(pd.DataFrame(ns),'results/nonrank_direction_summary.tsv');save(pd.DataFrame(bs),'results/random_reference_summary.tsv')

gsea=read(BASE/'results/all_ligand_programme_responses.tsv');comp=read(ROOT/'results/composition_standardized_expression.tsv');pool=read(BASE/'inputs/all_ligand_pool_eligibility.tsv');original=read(BASE/'results/response_pattern_summary.tsv').set_index(['ligand','celltype']);classrows=[]
for r in inv.itertuples():
 b=effect[(effect.ligand==r.ligand)&(effect.celltype==r.celltype)];po=pool[(pool.ligand==r.ligand)&(pool.celltype==r.celltype)]
 row=dict(ligand=r.ligand,family=r.celltype,eligible=r.eligibility=='ELIGIBLE',exclusion_reason='' if r.eligibility=='ELIGIBLE' else r.eligibility,condition_class=r.condition_class)
 donors={}
 for method,short in [('raw','raw'),('bg_relative','bg'),('nonrank','nonrank')]:
  w=b.pivot(index='donor',columns='module',values=method).reindex(index=DONORS,columns=MODS)
  row[method+'_class']=pattern(w[UP],w[DEP]);row[short+'_up_concordance']=int((w[UP]<0).sum()) if w[UP].notna().all() else np.nan;row[short+'_depleted_concordance']=int((w[DEP]>0).sum()) if w[DEP].notna().all() else np.nan
  donors[method]={d:{m:float(w.loc[d,m]) if pd.notna(w.loc[d,m]) else None for m in MODS} for d in DONORS}
  if r.eligibility!='ELIGIBLE':row[method+'_class']='NE'
 row['continuous_effects_by_donor']=json.dumps(donors,separators=(',',':'),allow_nan=False)
 row['original_both']=original.loc[(r.ligand,r.celltype),'pattern']=='BOTH'
 assert row['raw_class']==original.loc[(r.ligand,r.celltype),'pattern'],(r.ligand,r.celltype)
 row['raw_and_nonrank_both']=row['raw_class']==row['nonrank_class']=='BOTH';row['raw_bg_nonrank_both']=row['raw_class']==row['bg_relative_class']==row['nonrank_class']=='BOTH'
 row['new_bg_both']=row['bg_relative_class']=='BOTH' and not row['original_both'];row['new_nonrank_both']=row['nonrank_class']=='BOTH' and not row['original_both']
 row['control_ids']=';'.join(sorted({s for ss in b.control_ids.unique() for s in ss.split(';')}))
 row['pool_coverage']=json.dumps({d:sorted(po[(po.donor==d)&po.eligible].pool.unique()) for d in DONORS},separators=(',',':'))
 c=comp[(comp.ligand==r.ligand)&(comp.celltype==r.celltype)]
 row['composition_assessment_status']='OUTSIDE_PRESPECIFIED_CONTEXT_SUBSET' if r.ligand not in SELECTED else 'NO_ELIGIBLE_PAIRED_CELLS' if c.empty else 'NO_FINER_RELEASED_SUBSTATE' if (c.annotation_resolution=='NO_FINER_RELEASED_SUBSTATE').all() else 'ESTIMATED_ALL_DONORS' if c.donor.nunique()==3 and (c.status=='ESTIMATED').all() else 'PARTIAL_COMMON_SUPPORT'
 gs=gsea[(gsea.ligand==r.ligand)&(gsea.celltype==r.celltype)]
 for dest,col in [('original_GSEA_NES','NES'),('original_GSEA_q','FDR_atlas'),('GSEA_test_family','FDR_scope')]:
  row[dest]=json.dumps({v.module:None if pd.isna(getattr(v,col)) else getattr(v,col) for v in gs.itertuples()},separators=(',',':'),allow_nan=False)
 classrows.append(row)
cl=pd.DataFrame(classrows);main=cl[cl.condition_class=='ligand'].copy();save(main,'results/classification_stability.tsv');save(cl[cl.condition_class=='assay_reference'],'results/No_Cytomix_classification.tsv');save(main[main.original_both],'results/original_eight_Both_tracking.tsv');save(main[main.new_bg_both|main.new_nonrank_both],'results/additional_Both_combinations.tsv')
joint=[];marginal=[];counts=[];double=[]
for f in ['ALL']+FAMS:
 b=main if f=='ALL' else main[main.family==f];e=b[b.eligible]
 for u,d in itertools.product(range(4),repeat=2):joint.append(dict(family=f,up_down_donor_count=u,depleted_up_donor_count=d,n_combinations=int(((e.raw_up_concordance==u)&(e.raw_depleted_concordance==d)).sum())))
 for axis,col in [('up_down','raw_up_concordance'),('depleted_up','raw_depleted_concordance')]:
  for n in range(4):marginal.append(dict(family=f,axis=axis,donor_count=n,n_combinations=int((e[col]==n).sum())))
 for me in ['raw','bg_relative','nonrank']:
  for cat in ['BOTH','UP_ONLY','DEPLETED_ONLY','OTHER_OBSERVED','NE']:counts.append(dict(family=f,method=me,classification=cat,n_combinations=int((b[me+'_class']==cat).sum()),n_planned=len(b),n_eligible=int(b.eligible.sum())))
 eb=eligible if f=='ALL' else eligible[eligible.celltype==f]
 for donor in ['ALL']+DONORS:
  q=eb if donor=='ALL' else eb[eb.donor==donor];w=q.pivot(index=['ligand','celltype','donor'],columns='module',values='raw')
  double.append(dict(family=f,donor=donor,n_matched_pairs=len(w),same_donor_both_decline=int(((w[UP]<0)&(w[DEP]<0)).sum()),same_donor_anti_disease=int(((w[UP]<0)&(w[DEP]>0)).sum()),same_donor_disease_direction=int(((w[UP]>0)&(w[DEP]<0)).sum())))
save(pd.DataFrame(joint),'results/raw_joint_direction_table.tsv');save(pd.DataFrame(marginal),'results/raw_marginal_direction_counts.tsv');save(pd.DataFrame(counts),'results/classification_counts_by_method.tsv');save(pd.DataFrame(double),'results/raw_same_donor_joint_counts.tsv')
ref_expected={'up_down':[35,61,87,95],'depleted_up':[143,89,31,15]};marg=pd.DataFrame(marginal);version={}
for axis,expected in ref_expected.items():
 observed=marg[(marg.family=='ALL')&(marg.axis==axis)].sort_values('donor_count').n_combinations.tolist();version[axis]=dict(reference=expected,observed=observed,match=observed==expected)
version['all_up_down_no_dep_up']=int(((main.raw_up_concordance==3)&(main.raw_depleted_concordance==0)).sum());version['Both']=int(main.original_both.sum());jw(ROOT/'review/raw_version_comparison.json',version)

# Structural ligand group from the task, not a receptor/SMAD-pathway equivalence.
membership=[]
for ligand in sorted(main.ligand.unique()):
 sub='BMP' if ligand in TGF and ligand.startswith('BMP') else 'GDF' if ligand in TGF and ligand.startswith('GDF') else 'TGFB' if ligand in TGF and ligand.startswith('TGFB') else 'INHBA' if ligand=='INHBA' else 'Other'
 membership.append(dict(ligand=ligand,TGF_beta_superfamily=ligand in TGF,subgroup=sub,source=str(TASK)+' section E1',grouping='task-specified structural ligand superfamily; not shared receptor or SMAD pathway'))
mem=pd.DataFrame(membership);save(mem,'results/ligand_family_membership.tsv');main=main.merge(mem,on='ligand')
fr=[]
for fam in ['ALL']+FAMS:
 for gr in ['TGF_beta_superfamily','Other']:
  b=main[main.TGF_beta_superfamily==(gr=='TGF_beta_superfamily')];b=b if fam=='ALL' else b[b.family==fam]
  for method in ['raw','bg_relative','nonrank']:
   for cat in ['BOTH','UP_ONLY','DEPLETED_ONLY','OTHER_OBSERVED','NE']:
    hit=b[b[method+'_class']==cat]
    fr.append(dict(family=fam,ligand_group=gr,method=method,classification=cat,n_combinations=len(hit),n_planned_combinations=len(b),n_eligible_combinations=int(b.eligible.sum()),n_ligands_with_category=hit.ligand.nunique(),n_ligands_planned=b.ligand.nunique(),n_ligands_with_eligible_family=b.loc[b.eligible,'ligand'].nunique(),ligands=';'.join(sorted(hit.ligand.unique()))))
save(pd.DataFrame(fr),'results/family_response_summary.tsv')
refs=effect[effect.ligand.isin(['IFNG','TNFA'])].copy();refs['source_condition']=refs.ligand;refs['public_label']=refs.ligand.map({'IFNG':'IFN-gamma','TNFA':'TNF-alpha'});save(refs,'results/reference_perturbation_effects.tsv');save(main[main.ligand.isin(['IFNG','TNFA'])],'results/reference_perturbation_classification.tsv')
poolrows=read(ROOT/'results/pool_background_summary.tsv')[['ligand','celltype','donor','pool']].drop_duplicates();poolrows=poolrows.merge(mem,on='ligand',how='left');save(poolrows,'results/ligand_pool_family_design.tsv')
des=[]
for p,b in poolrows.groupby('pool'):
 des.append(dict(pool=p,n_ligands=b.ligand.nunique(),n_TGF_ligands=b.loc[b.TGF_beta_superfamily==True,'ligand'].nunique(),TGF_ligands=';'.join(sorted(b.loc[b.TGF_beta_superfamily==True,'ligand'].unique())),other_ligands=';'.join(sorted(b.loc[b.TGF_beta_superfamily==False,'ligand'].unique()))))
save(pd.DataFrame(des),'results/pool_family_design_summary.tsv')

ov=[]
for a,b in itertools.combinations(MODS,2):
 for scope,A,B in [('frozen_full',set(full[a]),set(full[b])),('measured',set(sets[a]),set(sets[b]))]:
  common=A&B;names=sorted(common) if scope=='frozen_full' else sorted(genes[int(i)] for i in common)
  ov.append(dict(programme_a=a,programme_b=b,scope=scope,n_a=len(A),n_b=len(B),n_intersection=len(common),Jaccard=len(common)/len(A|B),intersection_fraction_a=len(common)/len(A),intersection_fraction_b=len(common)/len(B),common_genes=';'.join(names)))
save(pd.DataFrame(ov),'results/gene_set_overlap.tsv')
exc=read(ROOT/'results/overlap_removed_donor_effects.tsv');pairs=read(ROOT/'inputs/overlap_sensitivity_pairs.tsv');cor=[]
for definition,data in [('fixed_members',effect),('symmetric_overlap_removed',exc)]:
 data=data[data.condition_class=='ligand']
 for pair in pairs.itertuples():
  aa=pair.left if definition=='fixed_members' else pair.left_exclusive;bb=pair.right if definition=='fixed_members' else pair.right_exclusive
  for fam in FAMS:
   for donor in ['MEAN_OF_3']+DONORS:
    b=data[data.celltype==fam];b=b.groupby(['ligand','module'],as_index=False)[['raw','bg_relative','nonrank']].mean() if donor=='MEAN_OF_3' else b[b.donor==donor]
    for method in ['raw','bg_relative','nonrank']:
     w=b[b.module.isin([aa,bb])].pivot(index='ligand',columns='module',values=method).reindex(columns=[aa,bb]).dropna();valid=len(w)>=3 and w[aa].nunique()>1 and w[bb].nunique()>1
     cor.append(dict(membership=definition,celltype=fam,donor=donor,method=method,programme_a=pair.left,programme_b=pair.right,n_common_ligands=len(w),common_ligands=';'.join(sorted(w.index)),Pearson=w[aa].corr(w[bb]) if valid else np.nan,Spearman=spearmanr(w[aa],w[bb]).statistic if valid else np.nan,status='DESCRIPTIVE_ESTIMATED' if valid else 'NE_INSUFFICIENT_COVERAGE_OR_VARIATION',n_a_remaining=pair.n_left_remaining if definition!='fixed_members' else len(sets[aa]),n_b_remaining=pair.n_right_remaining if definition!='fixed_members' else len(sets[bb]),P_value='NOT_COMPUTED_SHARED_DONORS_AND_CONTROLS'))
save(pd.DataFrame(cor),'results/reference_correlations_by_method.tsv')
print('SUMMARY_COMPLETE',len(main),'planned combinations',main[['raw_class','bg_relative_class','nonrank_class']].apply(pd.Series.value_counts).to_string(),flush=True)
