import os
from pathlib import Path
import json,hashlib,itertools,shutil
import numpy as np,pandas as pd
from scipy.stats import norm,rankdata
from scipy.io import mmread
O=Path(__file__).resolve().parents[1];E=O/'organoid_donor_extension'
files=[Path(os.environ['ICI_ORGANOID_COMPLETE'])];assert files[0].is_file()
complete=json.loads(files[-1].read_text());R=Path(complete['output']);plan=json.loads((E/'review/analysis_plan_locked.json').read_text())
ref=plan['assay_reference'].split(';')[0].strip();conditions=[plan['primary_ligand']]+plan['secondary_ligands']+[ref]
r=pd.read_csv(R/'module_GSEA.tsv',sep='\t');meta=pd.read_csv(E/'source_data/pseudobulk_metadata.tsv',sep='\t')
assert set(r.ligand)==set(conditions) and len(r)==len(conditions)*len(plan['celltypes'])*len(plan['modules'])
assert not r.duplicated(['ligand','celltype','module']).any()
for ligand,ct in itertools.product(conditions,plan['celltypes']):
 n=meta[(meta.ligand==ligand)&(meta.celltype==ct)].donor.nunique();z=r[(r.ligand==ligand)&(r.celltype==ct)]
 assert (z.n_donors==n).all()
 assert z.status.eq('ESTIMATED').all() if n>=3 else z.NES.isna().all()
def bh(p,n):
 p=np.asarray(p);q=np.full(len(p),np.nan);ix=np.flatnonzero(np.isfinite(p));order=ix[np.argsort(p[ix])];q[order]=np.minimum(np.minimum.accumulate((p[order]*n/np.arange(1,len(order)+1))[::-1])[::-1],1);return q
for family,z in r.groupby('family'):
 assert np.allclose(bh(z.pval,len(z)),z.FDR_predeclared_family,equal_nan=True,atol=1e-12)
model=pd.read_csv(R/'model_inventory.tsv',sep='\t');gene_checks=[]
for z in model.itertuples():
 f=R/'gene_results'/f'{z.ligand}__{z.celltype}.tsv.gz';de=pd.read_csv(f,sep='\t');valid=de.stat.notna()
 assert not de.gene.duplicated().any()
 assert np.allclose(de.loc[valid,'stat'],de.loc[valid,'log2FoldChange']/de.loc[valid,'lfcSE'],rtol=1e-11,atol=1e-11)
 assert np.allclose(de.loc[valid,'pvalue'],2*norm.sf(abs(de.loc[valid,'stat'])),rtol=1e-9,atol=1e-13)
 assert np.allclose(de.padj,bh(de.pvalue,int(de.pvalue.notna().sum())),equal_nan=True,atol=1e-12)
 gene_checks.append(dict(ligand=z.ligand,celltype=z.celltype,rows=len(de),finite_Wald=int(valid.sum()),status='PASS'))
scores=pd.read_csv(R/'donor_module_scores.tsv',sep='\t');pairs=pd.read_csv(R/'donor_paired_effects.tsv',sep='\t')
counts=mmread(E/'source_data/pseudobulk_counts.mtx').tocsc();genes=pd.read_csv(E/'source_data/pseudobulk_genes.tsv',sep='\t').gene
defs=pd.read_csv(O/'supplementary/original_tables/Table_S21_frozen_module_definitions.tsv',sep='\t')
geneindex={g:i for i,g in enumerate(genes)};scores_checked=0
for z in scores.itertuples():
 j=int(meta.index[meta.sample_id==z.sample_id][0]);v=counts[:,j].toarray().ravel();ranks=rankdata(v)/len(v);ix=[geneindex[g] for g in defs.loc[defs.module==z.module,'gene'].unique() if g in geneindex];expected=float(ranks[ix].mean())
 assert np.isclose(z.raw_mean_rank,expected,rtol=0,atol=1e-13);scores_checked+=1
for z in pairs.itertuples():
 s=scores[(scores.ligand==z.ligand)&(scores.celltype==z.celltype)&(scores.module==z.module)].pivot(index='donor',columns='treatment',values='raw_mean_rank');v=(s.Treatment-s.Control).to_numpy()
 null=[np.mean(v*np.array(signs)) for signs in itertools.product([-1,1],repeat=len(v))];p=np.mean(np.abs(null)>=abs(v.mean())-1e-15)
 assert np.isclose(p,z.exact_signflip_P,atol=1e-14) and np.isclose(v.mean(),z.mean_paired_rank_change,atol=1e-14)
for family,z in pairs.groupby('family'):
 n=int(r.family.eq(family).sum());assert np.allclose(bh(z.exact_signflip_P,n),z.FDR_predeclared_family,atol=1e-12)
# Compare repeated unaffected estimates with the first attempt; fixes cannot change them.
old=pd.read_csv(Path(os.environ['ICI_ORGANOID_COMPARISON_GSEA']),sep='\t');common=r[r.ligand!=ref].merge(old,on=['ligand','celltype','module'],suffixes=('_new','_first'),validate='one_to_one')
for c in ['NES','pval','FDR_predeclared_family']:assert np.allclose(common[c+'_new'],common[c+'_first'],equal_nan=True,atol=1e-12)
pd.DataFrame(gene_checks).to_csv(E/'review/gene_model_numerical_checks.tsv',sep='\t',index=False)
old_g=pd.read_csv(O/'source_data/panel_sources/F4B.tsv',sep='\t')
comparison=r[r.ligand=='TNFSF12'].merge(old_g[['celltype','module','NES','pval','screen_padj']],on=['celltype','module'],how='left',suffixes=('_donor','_V1'))
comparison.to_csv(E/'source_data/TNFSF12_V1_vs_donor_models.tsv',sep='\t',index=False)
core=pd.read_csv(O/'review/reciprocity_core_gene_evidence.tsv',sep='\t');coregenes=sorted(core.gene.unique());coreout=[];sorl=[]
for ct in plan['celltypes']:
 f=R/'gene_results'/f'TNFSF12__{ct}.tsv.gz';de=pd.read_csv(f,sep='\t');de['source_file']=str(f);de['source_row']=np.arange(1,len(de)+1)
 coreout.append(de[de.gene.isin(coregenes)].copy());sorl.append(de[de.gene=='SORL1'].copy())
coreout=pd.concat(coreout,ignore_index=True);sorl=pd.concat(sorl,ignore_index=True)
assert len(coreout)==len(coregenes)*len(plan['celltypes'])
coreout.to_csv(E/'source_data/frozen_core_donor_gene_effects.tsv',sep='\t',index=False);sorl.to_csv(E/'source_data/SORL1_donor_gene_effects.tsv',sep='\t',index=False)
summary=dict(status='PASS',result_directory=str(R),planned_tests=len(r),estimated_tests=int(r.status.eq('ESTIMATED').sum()),models=len(model),numeric_gene_rows=sum(x['finite_Wald'] for x in gene_checks),scores_checked=scores_checked,exact_tests_checked=len(pairs),unchanged_primary_secondary_tests=len(common),primary_GSEA_FDR05=int(((r.family=='primary_TNFSF12')&(r.FDR_predeclared_family<.05)).sum()),primary_donor_exact_FDR05=int(((pairs.family=='primary_TNFSF12')&(pairs.FDR_predeclared_family<.05)).sum()),core_TA_positive=int(((coreout.celltype=='Transit_Amplifying')&(coreout.log2FoldChange>0)).sum()),core_TA_gene_FDR05=int(((coreout.celltype=='Transit_Amplifying')&(coreout.padj<.05)).sum()))
(E/'review/numerical_verification.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
# Preserve the first attempt completion record before promoting the verified run.
canonical=E/'review/analysis_complete.json'
if canonical.exists():canonical.rename(E/'review/first_attempt_runtime_completion.json')
shutil.copy2(files[-1],canonical)
print(json.dumps(summary,indent=2));print(r[r.ligand==ref][['celltype','module','NES','FDR_predeclared_family']].to_string(index=False))
