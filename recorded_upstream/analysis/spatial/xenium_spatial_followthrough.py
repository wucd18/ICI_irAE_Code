import os
"""Spatial follow-through in overlapping V1 donors, never an independent validation."""
from pathlib import Path
import json,sys,hashlib,h5py,numpy as np,pandas as pd,scipy
from scipy import sparse,stats
from datetime import datetime,timezone
from remote_h5_metadata import Remote,array
O=Path(__file__).resolve().parents[1];X=O/'xenium_extension'
if '--output-root' in sys.argv:X=Path(sys.argv[sys.argv.index('--output-root')+1])
for d in ['review','source_data','figures']:(X/d).mkdir(parents=True,exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
def rankscore(a,indices):
 ranks=stats.rankdata(a,axis=0,method='average')/a.shape[0]
 return np.vstack([ranks[i].mean(axis=0) for i in indices])
def fixture():
 a=np.array([[0,4],[2,0],[2,4],[0,1]],float);idx=[np.array([1,2])]
 assert np.array_equal(rankscore(a,idx),rankscore(np.log2(a/a.sum(axis=0)*1e6+1),idx))
 z=sparse.csr_matrix(a.T);assert np.array_equal(np.asarray(z.sum(axis=0)).ravel(),a.sum(axis=1))
 print('PASS: tied ranks and deterministic count aggregation fixtures.',flush=True)
fixture()
if '--test' in sys.argv:raise SystemExit(0)
meta=O/'new_data/xenium5k_metadata'
obs=pd.read_csv(meta/'obs.tsv.gz',sep='\t',keep_default_na=False,low_memory=False)
var=pd.read_csv(meta/'var.tsv.gz',sep='\t',keep_default_na=False)
mapping={'240927_Condition':'condition','25_06_11_Patient_ID_HS':'patient','25_06_11_ICI_5K_Compartments':'compartment','24_11_12_ICI_5K_Fine_annotations':'celltype','_index':'cell_key'}
obs=obs.rename(columns=mapping);assert not obs.cell_key.duplicated().any() and not var._index.duplicated().any()
frozen=Path(os.environ['ICI_FROZEN_MODULES'])
defs=pd.read_csv(frozen,sep='\t')
mods=['CURATED_IFN_VISIBILITY','ICI_COLITIS_EPITHELIAL_UP_STRICT','ICI_COLITIS_EPITHELIAL_UP_NON_IFN','ICI_COLITIS_EPITHELIAL_LOSS_STRICT']
targets=['SORL1','TNFSF12','TNFRSF12A']
idx=[];coverage=[];gene_map=[]
for m in mods:
 g=defs.loc[defs.module.eq(m),'gene'].unique();matched=np.flatnonzero(var._index.isin(g));assert len(matched)>=10
 idx.append(matched);coverage.append({'module':m,'requested':len(g),'available':len(matched),'fraction':len(matched)/len(g),'panel_genes':len(var)})
 gene_map.extend({'module':m,'gene':v,'status':'PRESENT' if v in set(var._index) else 'NOT_IN_PANEL'} for v in g)
assert set(targets)<=set(var._index)
plan={'utc':datetime.now(timezone.utc).isoformat(),'scope':'spatial localization and patient-compartment descriptive estimates in overlapping V1 discovery donors; not independent validation','modules':mods,'targets':targets,'reason':'Only nominated original V1 marker and upstream ligand/receptor; selected before reading expression layer','gene_matching':'exact released gene symbols; all available members retained, no refitting or external imputation','gene_background':'complete released filtered targeted panel; dimensions derived from file','scores':'original mean within-pseudobulk rank; disease-depleted=1-raw; raw ranks for cell maps','inclusion':'all source-QC-passing cells with a assigned patient and condition; published compartment labels; no new clustering or manual ROI','replicate':'patient; all cells aggregated within patient and source compartment','statistics':'descriptive patient estimates and all patient points; no P values or new inferential tests because cohorts overlap; sample coverage comes from released metadata','limits':['partial module coverage on targeted panel','overlap with discovery donors must be checked against the released patient mapping','cell subtype and segmentation differences remain','no inference of ligand signaling or functional restoration from localization'],'coverage':coverage,'source_hashes':{str(p):sha(p) for p in [meta/'obs.tsv.gz',meta/'var.tsv.gz',frozen]}}
lock=X/'review/analysis_plan_locked.json'
if not lock.exists():lock.write_text(json.dumps(plan,indent=2))
else:
 prior=json.loads(lock.read_text());assert prior['source_hashes']==plan['source_hashes'] and prior['modules']==mods
pd.DataFrame(coverage).to_csv(X/'source_data/module_coverage.tsv',sep='\t',index=False)
pd.DataFrame(gene_map).to_csv(X/'source_data/frozen_gene_map.tsv',sep='\t',index=False)
var.to_csv(X/'source_data/panel_genes.tsv',sep='\t',index=False)
summary=obs.groupby(['patient','condition','compartment'],dropna=False).size().reset_index(name='source_cells')
summary.to_csv(X/'source_data/patient_compartment_manifest.tsv',sep='\t',index=False)
if '--preflight' in sys.argv:
 print(json.dumps(plan,indent=2),flush=True);raise SystemExit(0)
source=next(f for f in json.loads((O/'new_data/figshare27327813.json').read_text())['files'] if f['name']=='25_11_12_Xenium_Dataset2_5K_Annotated.h5ad')
matrixfile=O/'new_data/xenium5k_metadata/released_raw_counts.npz'
if not matrixfile.exists():
 rem=Remote(source['download_url'],source['size'],meta/'byte_ranges')
 with h5py.File(rem,'r') as h:
  g=h['layers/raw_counts'];enc=g.attrs['encoding-type'];print('RAW layer',enc,dict(g.attrs),flush=True)
  assert enc in ['csr_matrix','csc_matrix']
  arrays={}
  for k in ['data','indices','indptr']:
   ds=g[k];print(k,ds.shape,ds.dtype,flush=True)
   arrays[k]=np.concatenate([ds[j:j+10_000_000] for j in range(0,len(ds),10_000_000)])
  ctor=sparse.csr_matrix if enc=='csr_matrix' else sparse.csc_matrix
  a=ctor((arrays['data'],arrays['indices'],arrays['indptr']),shape=tuple(g.attrs['shape'])).tocsr()
 sparse.save_npz(matrixfile,a)
 (meta/'counts_acquisition.json').write_text(json.dumps({'url':source['download_url'],'remote_bytes':source['size'],'published_remote_md5':source.get('computed_md5'),'local_raw_count_npz_sha256':sha(matrixfile),'downloaded_ranges':{str(v['start']):v for v in rem.records},'utc':datetime.now(timezone.utc).isoformat(),'scope':'released layers/raw_counts only, preserves sparse values and order; not whole annotated object'},indent=2))
else:a=sparse.load_npz(matrixfile)
assert a.shape==(len(obs),len(var)),(a.shape,len(obs),len(var))
assert np.isfinite(a.data).all() and (a.data>=0).all() and np.equal(a.data,np.floor(a.data)).all()
lib=np.asarray(a.sum(axis=1)).ravel()
# Source obs QC values can precede gene filtering; never fill removed molecules/genes.
qc_delta=obs.total_counts.to_numpy()-lib
assert np.isfinite(qc_delta).all() and (qc_delta>=0).all() and np.equal(qc_delta,np.floor(qc_delta)).all(),'Released panel has counts not present in source QC totals'
audit=obs.loc[qc_delta!=0,['cell_key','patient','condition','total_counts']].copy()
audit['released_panel_count_sum']=lib[qc_delta!=0];audit['QC_minus_panel']=qc_delta[qc_delta!=0]
audit.to_csv(X/'review/source_QC_vs_filtered_panel.tsv',sep='\t',index=False)
(X/'review/source_QC_resolution.json').write_text(json.dumps({'source':'https://www.jci.org/articles/view/202488; Data preprocessing and annotation','source_rule':'Genes with <1 count and detected in <10 cells removed in Datasets 1 and 2','finding':'Source obs QC totals exceed the released filtered-panel count sum in the listed cells; compatible with source gene filtering, exact omitted gene identity not inferred','action':'Use only the released raw_counts layer, compute denominators from that layer; no missing count added; all source-accepted cells retained','mismatch_cells':len(audit),'total_difference':float(qc_delta.sum()),'statistical_inputs_changed':False},indent=2))
use=(~obs.patient.eq('unassigned'))&(~obs.condition.eq('unassigned'))&(lib>0)
obs['analysis_included']=use
cols=['cell_key','cell_id','x_centroid','y_centroid','patient','condition','compartment','celltype','analysis_included']
celltable=obs[cols].copy();out=[];genes=[];pb=[];pbkeys=[]
for (patient,condition,comp),rows in obs[use].groupby(['patient','condition','compartment']).groups.items():
 b=a[np.asarray(rows),:];v=np.asarray(b.sum(axis=0)).ravel();s=rankscore(v[:,None],idx)[:,0];key=f'{patient}|{comp}'
 pb.append(v);pbkeys.append({'record_key':key,'patient':patient,'condition':condition,'compartment':comp,'n_cells':len(rows)})
 for i,m in enumerate(mods):
  out.append({'patient':patient,'condition':condition,'compartment':comp,'n_cells':len(rows),'module':m,'raw_rank':s[i],'oriented_rank':1-s[i] if m.endswith('LOSS_STRICT') else s[i],'available_genes':len(idx[i]),'record_key':key+'|'+m})
 for gene in targets:
  gi=int(np.flatnonzero(var._index.eq(gene))[0]);z=b[:,gi].toarray().ravel()
  genes.append({'patient':patient,'condition':condition,'compartment':comp,'gene':gene,'n_cells':len(rows),'total_molecules':float(z.sum()),'mean_molecules_per_cell':float(z.mean()),'fraction_positive_cells':float((z>0).mean()),'CPM':float(z.sum()/v.sum()*1e6),'record_key':key+'|'+gene})
for gene in targets:
 gi=int(np.flatnonzero(var._index.eq(gene))[0]);celltable[gene+'_molecules']=a[:,gi].toarray().ravel()
for j in range(0,len(obs),1024):
 z=rankscore(a[j:j+1024,:].toarray().T,idx)
 for i,m in enumerate(mods):celltable.loc[j:j+z.shape[1]-1,m]=z[i]
celltable.to_csv(X/'source_data/cell_coordinates_scores.tsv.gz',sep='\t',index=False,compression='gzip')
pd.DataFrame(out).to_csv(X/'source_data/patient_compartment_scores.tsv',sep='\t',index=False)
pd.DataFrame(genes).to_csv(X/'source_data/patient_gene_expression.tsv',sep='\t',index=False)
pd.DataFrame(pbkeys).to_csv(X/'source_data/pseudobulk_column_keys.tsv',sep='\t',index=False)
sparse.save_npz(X/'source_data/patient_compartment_counts.npz',sparse.csc_matrix(np.column_stack(pb)))
result={'utc':datetime.now(timezone.utc).isoformat(),'status':'EXECUTED_DESCRIPTIVE','included_cells':int(use.sum()),'excluded_cells':int((~use).sum()),'patients':int(obs[use].patient.nunique()),'patient_counts':obs[use].drop_duplicates('patient').condition.value_counts().to_dict(),'source_raw_npz_sha256':sha(matrixfile),'hypothesis_tests':0,'independent_new_cohorts':0,'software':{'python':sys.version,'scipy':scipy.__version__,'h5py':h5py.__version__},'source_hashes_verified':all(sha(p)==v for p,v in plan['source_hashes'].items())}
assert result['source_hashes_verified'];(X/'review/analysis_summary.json').write_text(json.dumps(result,indent=2));print(result,flush=True)
