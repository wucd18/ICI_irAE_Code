from pathlib import Path
import json,hashlib,sys,itertools
import numpy as np,pandas as pd
from scipy.sparse import csc_matrix
from scipy.io import mmwrite
from datetime import datetime,timezone
from stream_r_object import Node
from export_r_capture import unpack
sys.stdout.reconfigure(encoding='utf-8')
O=Path(__file__).resolve().parents[1];D=O/'new_data/organoid_extracted';E=O/'organoid_donor_extension'
for n in ['source_data','review','results','figures']:(E/n).mkdir(parents=True,exist_ok=True)
arrays,obs,genes,shape=unpack(D)
nom=pd.read_csv(O/'source_data/panel_sources/F4A.tsv',sep='\t')
ligands=sorted(nom.loc[nom['class'].eq('Directional nominees'),'perturbation'].tolist())
assert len(ligands)==len(set(ligands)) and 'TNFSF12' in ligands
conditions=ligands+['No_Cytomix'];cts=['Stem','Transit_Amplifying','Colonocyte','Goblet'];controls=['Sato_1','Sato_2','Sato_3']
mapping={'Colonocyte_Immature':'Colonocyte','Colonocyte_Mature':'Colonocyte','Goblet_Immature':'Goblet','Goblet_Mature':'Goblet'}
assert obs.final_annotations.replace(mapping).equals(obs.celltype)
assert obs[['assignment','perturbation','pool','celltype']].notna().all().all()
assert set(obs.demux_type)=={'singlet'}
download=json.loads((O/'new_data/organoid_download/download_complete.json').read_text(encoding='utf-8'))
plan=dict(utc=datetime.now(timezone.utc).isoformat(),source=download,source_count_slot='assays/RNA/counts; complete serialized gene background',celltypes=cts,source_annotation='Released celltype; confirmed identical to final_annotations with mature/immature colonocytes and goblets merged as source Figure 4 code lines 333-337.',primary_ligand='TNFSF12',secondary_ligands=[x for x in ligands if x!='TNFSF12'],assay_reference='No_Cytomix; not a newly screened candidate',control=controls,matching='Same donor, pool and published celltype; pool retained only if treatment and combined cytomix-only controls each contain >=20 released cells. Sum counts across eligible matched pools and technical control wells within donor before inference.',unit='donor; require >=3 matched donors per condition-celltype, no cell or pool pseudoreplication',modules=['ICI_COLITIS_EPITHELIAL_UP_STRICT','ICI_COLITIS_EPITHELIAL_UP_NON_IFN','ICI_COLITIS_EPITHELIAL_LOSS_STRICT','CURATED_IFN_VISIBILITY'],gene_prefilter='Condition-blind total raw count across paired pseudobulks >=10; no effect or P-value filter.',DE='DESeq2 1.44.0 design ~ donor + treatment; no LFC shrinkage, no automatic outlier replacement; independentFiltering=FALSE; retain Cook outlier flags and exclude nonfinite/P-unavailable ranks.',enrichment='Unshrunken Wald statistic; fgseaSimple, 10000 permutations, minSize=10, maxSize=500, scoreType=std, nproc=1, seed=20260826. All finite ranks, no significance prefilter.',multiple_testing='BH: all 16 TNFSF12 celltype-module tests as primary; all 96 other original-nominee tests as secondary; 16 No_Cytomix assay-reference tests separately. Planned missing estimates remain NA and retain family denominator.',donor_scores='Mean within-donor/condition gene rank using all released genes; positive raw-rank change in disease-depleted module is favourable. Exact two-sided sign-flip P from paired donor changes, with same predeclared multiplicity families. This has low power with three donors.',scope='Targeted reanalysis of original nominations, not a new candidate screen or independent perturbation replication. Preserve complete original estimates. No original threshold, set or candidate rank is overwritten.',limitations='Three healthy organoid donors, source injury model, no clinical safety or functional-rescue endpoint. GSEA tests gene-set enrichment and does not add donor replicates.')
with (E/'review/analysis_plan_locked.json').open('x',encoding='utf-8') as f:json.dump(plan,f,indent=2)
idx=obs.reset_index();idx['column_index']=np.arange(len(idx));idx['arm']=np.where(idx.perturbation.isin(controls),'Control',idx.perturbation)
g=idx.groupby(['assignment','pool','celltype','arm'],sort=True).size().rename('n_cells').reset_index()
g.to_csv(E/'source_data/all_released_group_cell_counts.tsv',sep='\t',index=False)
include=[];aud=[]
for ligand,ct,donor in itertools.product(conditions,cts,sorted(idx.assignment.unique())):
 eligible=[]
 for pool in sorted(idx.pool.unique()):
  block=g[(g.assignment==donor)&(g.pool==pool)&(g.celltype==ct)]
  ntr=int(block.loc[block.arm==ligand,'n_cells'].sum());nco=int(block.loc[block.arm=='Control','n_cells'].sum())
  if ntr:
   ok=ntr>=20 and nco>=20;aud.append(dict(ligand=ligand,celltype=ct,donor=donor,pool=pool,treatment_cells=ntr,control_cells=nco,eligible=ok,reason='PASS' if ok else 'one_arm_below_original_20_cell_minimum'))
   if ok:eligible.append(pool)
 if eligible:
  for arm in ['Control',ligand]:
   select=idx.assignment.eq(donor)&idx.pool.isin(eligible)&idx.celltype.eq(ct)&idx.arm.eq(arm)
   include.append(dict(ligand=ligand,celltype=ct,donor=donor,treatment='Control' if arm=='Control' else 'Treatment',pools=';'.join(eligible),n_cells=int(select.sum()),indices=np.flatnonzero(select)))
audit=pd.DataFrame(aud);audit.to_csv(E/'source_data/matched_pool_eligibility.tsv',sep='\t',index=False)
meta=pd.DataFrame([{k:v for k,v in z.items() if k!='indices'} for z in include]);meta.insert(0,'sample_id',[f'ORG_PB{i+1:04d}' for i in range(len(meta))])
meta.to_csv(E/'source_data/pseudobulk_metadata.tsv',sep='\t',index=False)
# Audit every released cell, independently against metadata, using bounded chunks.
p=np.asarray(arrays['p'],dtype=np.int64);x=arrays['x'];ii=arrays['i'];sums=np.zeros(shape[1]);features=np.diff(p)
for lo in range(0,shape[1],2048):
 hi=min(shape[1],lo+2048);start,end=p[lo],p[hi];v=np.asarray(x[start:end],dtype=np.float64);ri=np.asarray(ii[start:end],dtype=np.int64)
 assert np.isfinite(v).all() and (v>0).all() and np.equal(v,np.floor(v)).all()
 assert (ri>=0).all() and (ri<shape[0]).all()
 lengths=np.diff(p[lo:hi+1]);assert (lengths>0).all()
 sums[lo:hi]=np.add.reduceat(v,p[lo:hi]-start)
assert np.array_equal(sums,obs.nCount_RNA.to_numpy()),'Cell library sum mismatch'
assert np.array_equal(features,obs.nFeature_RNA.to_numpy()),'Cell feature count mismatch'
pd.DataFrame(dict(cell_id=obs.index,nCount_from_counts=sums,nCount_released=obs.nCount_RNA,nFeature_from_counts=features,nFeature_released=obs.nFeature_RNA)).to_csv(E/'source_data/cell_count_QC.tsv.gz',sep='\t',index=False,compression={'method':'gzip','compresslevel':1})
mat=np.zeros((shape[0],len(include)),dtype=np.int64);keys=[]
for j,z in enumerate(include):
 for cell in z['indices']:
  start,end=p[cell],p[cell+1];ri=np.asarray(ii[start:end],dtype=np.int64);v=np.asarray(x[start:end],dtype=np.int64)
  assert len(np.unique(ri))==len(ri),'Duplicate sparse row within cell'
  mat[:,j]+=np.bincount(ri,weights=v,minlength=shape[0]).astype(np.int64)
 keys.extend(dict(sample_id=meta.sample_id[j],cell_id=obs.index[cell],source_column_index=int(cell),ligand=z['ligand'],donor=z['donor'],pool=obs.pool.iloc[cell],celltype=z['celltype'],treatment=z['treatment']) for cell in z['indices'])
 assert mat[:,j].sum()==sums[z['indices']].sum()
 print('AGGREGATED',j+1,len(include),flush=True)
mmwrite(str(E/'source_data/pseudobulk_counts.mtx'),csc_matrix(mat),field='integer')
pd.DataFrame({'gene':genes}).to_csv(E/'source_data/pseudobulk_genes.tsv',sep='\t',index=False)
pd.DataFrame(keys).to_csv(E/'source_data/pseudobulk_cell_keys.tsv.gz',sep='\t',index=False,compression={'method':'gzip','compresslevel':1})
summary=dict(status='PASS',released_cells=shape[1],released_genes=shape[0],raw_nonzero_values=len(x),all_cell_nCount_exact=True,all_cell_nFeature_exact=True,pseudobulks=len(meta),unique_selected_cells=len(set(k['cell_id'] for k in keys)),source_serialized_object_sha256=download['sha256'])
(E/'review/preflight_QC.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(summary)
