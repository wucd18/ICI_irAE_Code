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

summary=dict(status="PASS",models=len(model),gene_rows_checked=sum(z["finite_Wald"] for z in gene_checks),donor_scores=scores_checked,exact_tests=len(pairs))
(E/"review/read_only_model_verification.json").write_text(json.dumps(summary,indent=2))
print(json.dumps(summary))
