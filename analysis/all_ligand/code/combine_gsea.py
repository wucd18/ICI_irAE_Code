from common import *
from analysis_functions import *
import itertools
cfg=jr(ROOT/'analysis_config.json');inv=pd.read_csv(ROOT/'inputs/comparison_inventory.tsv',sep='\t');jobs=pd.read_csv(ROOT/'inputs/new_model_jobs.tsv',sep='\t')
new=[];cov=[];models=[]
for r in jobs.itertuples():
    stem=ROOT/'results/new_models'/r.stem
    assert jr(str(stem)+'_complete.json')['status']=='COMPLETE',r.stem
    b=pd.read_csv(str(stem)+'_GSEA.tsv',sep='\t',float_precision='round_trip');b['value_origin']='FIT_NEW';b['result_source']=str(stem)+'_GSEA.tsv';new.append(b)
    cov.append(pd.read_csv(str(stem)+'_coverage.tsv',sep='\t'));models.append(pd.read_csv(str(stem)+'_model.tsv',sep='\t'))
old=pd.read_csv(ROOT/'inputs/targeted_module_GSEA.tsv',sep='\t',float_precision='round_trip')
old=old[old.status=='ESTIMATED'].rename(columns={'FDR_predeclared_family':'old_targeted_FDR','family':'old_targeted_FDR_family'});old['value_origin']='REUSED_TARGETED';old['result_source']=str(SAVED/'module_GSEA.tsv')
allg=pd.concat([old]+new,ignore_index=True);assert not allg.duplicated(['ligand','celltype','module']).any()
planned=pd.MultiIndex.from_product([cfg['ligands']+['No_Cytomix'],cfg['families'],cfg['modules_in_original_fgsea_order']],names=['ligand','celltype','module']).to_frame(index=False)
allg=planned.merge(allg.drop(columns=['n_donors']),how='left',validate='one_to_one').merge(inv[['ligand','celltype','n_donors','eligibility','condition_class','action']],validate='many_to_one')
allg['status']=allg.status.fillna(allg.eligibility);allg['value_origin']=allg.value_origin.fillna('NO_MODEL');allg['result_source']=allg.result_source.fillna(str(ROOT/'inputs/comparison_inventory.tsv'))
allg['FDR_scope']=np.where(allg.condition_class.eq('assay_reference'),'assay_reference_16',np.where(allg.module.eq(RES),'all_ligand_residual_316','all_ligand_primary_948'))
allg['FDR_planned_n']=allg.FDR_scope.map({'assay_reference_16':16,'all_ligand_residual_316':316,'all_ligand_primary_948':948})
allg['FDR_atlas']=np.nan
for scope,b in allg.groupby('FDR_scope'):
    n=int(b.FDR_planned_n.iloc[0]);assert len(b)==n;allg.loc[b.index,'FDR_atlas']=bh_planned(b.pval,n)
coverage=pd.concat([pd.read_csv(ROOT/'inputs/targeted_module_coverage.tsv',sep='\t')]+cov,ignore_index=True)
assert not coverage.duplicated(['ligand','celltype','module']).any()
allg=allg.merge(coverage,how='left',validate='one_to_one')
defs=pd.read_csv(ROOT/'inputs/frozen_module_definitions.csv');gm=set(pd.read_csv(ROOT/'inputs/pseudobulk_genes.tsv',sep='\t').gene)
for m,b in defs.groupby('module'):
    use=allg.module.eq(m);allg.loc[use,'requested']=allg.loc[use,'requested'].fillna(b.gene.nunique());allg.loc[use,'raw_available']=allg.loc[use,'raw_available'].fillna(b[b.gene.isin(gm)].gene.nunique())
stats=[]
delta=pd.read_csv(ROOT/'results/all_readout_donor_deltas_long.tsv',sep='\t',float_precision='round_trip')
for (l,c,m),b in delta.groupby(['ligand','celltype','module']):
    stats.append(dict(ligand=l,celltype=c,module=m,n_observed_score_pairs=len(b),delta_mean=b.delta.mean(),delta_median=b.delta.median(),delta_min=b.delta.min(),delta_max=b.delta.max(),n_delta_positive=int((b.delta>0).sum()),n_delta_negative=int((b.delta<0).sum()),n_delta_zero=int((b.delta==0).sum())))
allg=allg.merge(pd.DataFrame(stats),how='left',validate='one_to_one').merge(pd.read_csv(ROOT/'inputs/condition_metadata.tsv',sep='\t')[['ligand','dose','published_protein_name','mapping_status']],how='left',validate='many_to_one')
allg['raw_score_gene_background']=len(gm);allg['source_accession']='GSE313368';allg['score_sign']='Treatment minus matched control; disease-depleted not inverted'
allg.loc[allg.ligand=='IL22','mapping_status']='BLOCKED_DOSE_IL11_IL22_SOURCE_CONFLICT';allg.loc[allg.ligand=='No_Cytomix','mapping_status']='ASSAY_REFERENCE_NOT_A_LIGAND'
save(allg,'results/all_ligand_programme_responses.tsv');save(pd.concat(models),'results/new_model_inventory.tsv')
scopes=allg.groupby('FDR_scope').agg(planned=('module','size'),estimable=('pval','count'),significant=('FDR_atlas',lambda x:int((x<.05).sum()))).reset_index();scopes['NE']=scopes.planned-scopes.estimable
save(scopes,'results/FDR_scope_counts.tsv')
# Verify reused estimates by numerical identity, not a refit of old results.
reused=allg[allg.value_origin=='REUSED_TARGETED'].merge(old,on=['ligand','celltype','module'],suffixes=('_now','_old'),validate='one_to_one')
for col in ['NES','pval','ES','size','n_ranked_genes','old_targeted_FDR']:
    assert np.array_equal(reused[col+'_now'].to_numpy(),reused[col+'_old'].to_numpy(),equal_nan=True),col
assert len(reused)==100 and len(allg)==1280
assert allg.loc[allg.status=='ESTIMATED',['NES','pval','FDR_atlas']].notna().all().all()
assert allg.loc[allg.eligibility!='ELIGIBLE','FDR_atlas'].isna().all()
ref=allg[allg.condition_class=='assay_reference'];assert np.allclose(ref.FDR_atlas,ref.old_targeted_FDR,atol=1e-14,rtol=1e-12)
jw(ROOT/'review/result_checks.json',dict(status='PASS',models_new=len(jobs),models_reused=25,reused_module_estimates_unchanged=100,primary_planned=948,primary_estimable=int(scopes.loc[scopes.FDR_scope=='all_ligand_primary_948','estimable'].iloc[0]),planned_NE_retained=True,reference_FDR_preserved=True,generated_utc=now()))
print(scopes.to_string(index=False),flush=True)
