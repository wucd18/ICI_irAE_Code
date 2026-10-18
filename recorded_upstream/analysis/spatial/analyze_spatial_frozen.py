"""Bounded new-cohort transfer of fixed V1 scores; no discovery, GSEA or refitting."""
from pathlib import Path
import json,hashlib,sys,gzip
import numpy as np,pandas as pd,h5py,scipy
from scipy import sparse,stats,io
from datetime import datetime,timezone
O=Path(__file__).resolve().parents[1];X=O/'spatial_extension'
INPUT=X
if '--output-root' in sys.argv:X=Path(sys.argv[sys.argv.index('--output-root')+1])
for d in ['review','source_data','figures']:(X/d).mkdir(parents=True,exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def score(counts,indices):
 # log2(CPM+1) preserves within-column ordering and ties exactly.
 ranks=stats.rankdata(counts,axis=0,method='average')/counts.shape[0]
 return np.vstack([ranks[idx].mean(axis=0) for idx in indices])
def fixtures():
 a=np.array([[0,5],[0,0],[4,5],[2,1],[4,0]],float)
 i=[np.array([0,2,4]),np.array([1,3])]
 assert np.allclose(score(a,i),score(np.log2(a/a.sum(axis=0)*1e6+1),i),rtol=0,atol=0)
 assert score(a,i).shape==(2,2)
 assert np.allclose((1-score(a,i))+score(a,i),1)
 print('Fixtures PASS: tied zero ranks, monotone logCPM equivalence and orientation.',flush=True)
fixtures()
if '--test' in sys.argv:sys.exit(0)
plan=json.loads((INPUT/'review/analysis_plan_locked.json').read_text())
before=pd.read_csv(INPUT/'review/input_manifest_before.tsv',sep='\t');assert all(sha(r.path)==r.sha256 for r in before.itertuples())
meta=pd.read_csv(INPUT/'source_data/sample_manifest.tsv',sep='\t').fillna('')
defs=pd.read_csv(INPUT/'source_data/frozen_module_members.tsv',sep='\t')
features={r.section:pd.read_csv(INPUT/'source_data'/f'{r.section}_feature_index.tsv',sep='\t',keep_default_na=False) for r in meta.itertuples()}
if X!=INPUT:
 import shutil
 for rel in ['review/analysis_plan_locked.json','review/input_manifest_before.tsv','source_data/sample_manifest.tsv','source_data/frozen_module_members.tsv']:shutil.copy2(INPUT/rel,X/rel)
first=features[meta.iloc[0].section];common=set(first.feature_id)
for f in features.values():common&=set(f.feature_id)
universe=first[first.feature_id.isin(common)].copy();ids=universe.feature_id.to_numpy();symbols=universe.symbol.astype(str).to_numpy()
aliases={};alias_records=[]
for section,f in features.items():
 ff=f.set_index('feature_id').loc[ids]
 for gid,sym in zip(ids,ff.symbol.astype(str)):
  aliases.setdefault(sym,set()).add(gid)
  alias_records.append({'section':section,'feature_id':gid,'symbol':sym})
alias_table=pd.DataFrame(alias_records).drop_duplicates(['feature_id','symbol'])
alias_table.to_csv(X/'source_data/released_gene_aliases.tsv',sep='\t',index=False)
universe.to_csv(X/'source_data/common_feature_universe.tsv',sep='\t',index=False)
id_to_index={g:i for i,g in enumerate(ids)}
indices=[];coverage=[];mapping=[]
for m in plan['modules']:
 req=defs.loc[defs.module.eq(m),'gene'].astype(str).unique();idx=[]
 for gene in req:
  matches=aliases.get(gene,set())
  if len(matches)>1:raise ValueError(f'Ambiguous module gene {gene}: {matches}')
  gid=next(iter(matches)) if matches else ''
  mapping.append({'module':m,'original_gene':gene,'feature_id':gid,'status':'MAPPED' if gid else 'NOT_AVAILABLE'})
  if gid:idx.append(id_to_index[gid])
 assert len(idx)==len(set(idx)),('Duplicate gene identity within module',m)
 assert len(idx)>=10,m
 indices.append(np.asarray(idx));coverage.append(dict(module=m,requested=len(req),available=len(idx),missing=';'.join(g for g in req if g not in aliases),background_genes=len(ids)))
pd.DataFrame(mapping).to_csv(X/'source_data/frozen_gene_identity_map.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(X/'source_data/module_coverage.tsv',sep='\t',index=False)
expected={m:defs.loc[defs.module.eq(m),'expected_direction_in_irAE'].iloc[0] for m in plan['modules']}
patient_counts={};patient_spots={};spot_tables=[];qcs=[]
for r in meta.itertuples():
 if str(r.matrix).endswith('.h5'):
  with h5py.File(r.matrix) as f:
   g=f['matrix'];a=sparse.csc_matrix((g['data'][:],g['indices'][:],g['indptr'][:]),shape=tuple(g['shape'][:]));barcodes=[s.decode() for s in g['barcodes'][:]]
 else:
  with gzip.open(r.matrix,'rb') as f:a=io.mmread(f).tocsc()
  barcodes=pd.read_csv(r.barcodes,header=None,sep='\t')[0].tolist()
 assert len(set(barcodes))==len(barcodes)
 assert np.isfinite(a.data).all() and (a.data>=0).all() and np.equal(a.data,np.floor(a.data)).all()
 f=features[r.section];order=f.set_index('feature_id').loc[ids,'feature_row'].to_numpy()-1;a=a[order,:]
 pos=pd.read_csv(r.coordinates,header=None,names=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col']).set_index('barcode')
 assert not pos.index.duplicated().any() and set(barcodes)<=set(pos.index)
 pos=pos.loc[barcodes].copy();lib=np.asarray(a.sum(axis=0)).ravel();keep=(pos.in_tissue.to_numpy()==1)&(lib>0)
 a=a[:,keep];pos=pos.loc[np.asarray(barcodes)[keep]];lib=lib[keep]
 assert a.shape[1]>=20,r.section
 counts=np.asarray(a.sum(axis=1)).ravel().astype(np.int64)
 patient_counts[r.patient]=patient_counts.get(r.patient,np.zeros(len(ids),dtype=np.int64))+counts
 patient_spots[r.patient]=patient_spots.get(r.patient,0)+a.shape[1]
 raw=np.hstack([score(a[:,j:j+128].toarray(),indices) for j in range(0,a.shape[1],128)])
 scale=json.loads(Path(r.scalefactors).read_text());pos['x']=pos.pixel_col*scale['tissue_lowres_scalef'];pos['y']=pos.pixel_row*scale['tissue_lowres_scalef'];pos['library_count']=lib
 for i,m in enumerate(plan['modules']):
  d=pos.reset_index();d['section']=r.section;d['gsm']=r.gsm;d['patient']=r.patient;d['condition']=r.condition;d['module']=m;d['raw_score']=raw[i];d['oriented_score']=raw[i] if expected[m]=='up' else 1-raw[i];d['matrix_source']=r.matrix;d['record_key']=r.gsm+'|'+d.barcode
  spot_tables.append(d)
 qcs.append(dict(section=r.section,gsm=r.gsm,patient=r.patient,condition=r.condition,input_spots=len(barcodes),kept_tissue_spots=a.shape[1],zero_library_spots=int((np.asarray(a.sum(axis=0)).ravel()==0).sum()),total_UMI=int(lib.sum()),median_UMI=float(np.median(lib)),common_genes=len(ids),image=r.image,spot_diameter_lowres=scale['spot_diameter_fullres']*scale['tissue_lowres_scalef']))
 print(r.section,r.patient,a.shape,'completed',flush=True)
pd.concat(spot_tables,ignore_index=True).to_csv(X/'source_data/spot_scores.tsv.gz',sep='\t',index=False,compression='gzip')
qcs=pd.DataFrame(qcs);qcs.to_csv(X/'source_data/section_QC.tsv',sep='\t',index=False)
patients=list(patient_counts);mat=np.column_stack([patient_counts[p] for p in patients]);scores=score(mat,indices)
sparse.save_npz(X/'source_data/patient_pseudobulk_counts.npz',sparse.csc_matrix(mat))
pd.DataFrame({'patient':patients,'matrix_column_1based':np.arange(len(patients))+1}).to_csv(X/'source_data/patient_matrix_columns.tsv',sep='\t',index=False)
out=[]
for j,p in enumerate(patients):
 pm=meta[meta.patient.eq(p)];assert pm.condition.nunique()==1
 for i,m in enumerate(plan['modules']):
  raw=float(scores[i,j]);out.append(dict(patient=p,condition=pm.condition.iloc[0],sex=pm.sex.iloc[0],age=pm.age.iloc[0],sections=';'.join(pm.section),n_sections=len(pm),n_spots=patient_spots[p],module=m,expected_direction=expected[m],raw_score=raw,oriented_score=raw if expected[m]=='up' else 1-raw,available_genes=len(indices[i]),source_excel_row=int(pm.source_excel_row.iloc[0]),record_key=p+'|'+m))
out=pd.DataFrame(out);out.to_csv(X/'source_data/patient_scores.tsv',sep='\t',index=False)
tests=[]
for m,b in out.groupby('module',sort=False):
 x=b.loc[b.condition.eq('ICI colitis'),'oriented_score'].to_numpy();y=b.loc[b.condition.eq('Healthy'),'oriented_score'].to_numpy();assert min(len(x),len(y))>=plan['minimum_patients_per_group']
 t=stats.mannwhitneyu(x,y,alternative='two-sided');delta=(np.greater.outer(x,y).sum()-np.less.outer(x,y).sum())/(len(x)*len(y))
 tests.append(dict(module=m,n_case=len(x),n_control=len(y),median_case=float(np.median(x)),median_control=float(np.median(y)),median_difference=float(np.median(x)-np.median(y)),cliffs_delta=float(delta),U=float(t.statistic),P=float(t.pvalue),test='two-sided Mann-Whitney U, scipy default auto; exact when untied',record_key='Oxford|ICI_vs_Healthy|'+m))
tests=pd.DataFrame(tests);p=tests.P.to_numpy();order=np.argsort(p);q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];back=np.empty(len(p));back[order]=np.minimum(q,1);tests['FDR_new_four_module_family']=back
tests.to_csv(X/'source_data/frozen_transfer_results.tsv',sep='\t',index=False)
post=before.copy();post['sha256_after']=[sha(p) for p in post.path];post['status']=np.where(post.sha256==post.sha256_after,'PASS','FAIL');post.to_csv(X/'review/input_hash_comparison.tsv',sep='\t',index=False);assert post.status.eq('PASS').all()
summary={'utc':datetime.now(timezone.utc).isoformat(),'status':'EXECUTED','plan_sha256':sha(X/'review/analysis_plan_locked.json'),'sections':len(meta),'patients':len(patients),'tissue_spots':int(qcs.kept_tissue_spots.sum()),'common_gene_universe':len(ids),'duplicate_symbols':len(symbols)-len(set(symbols)),'module_coverage':coverage,'tests':tests.to_dict('records'),'software':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'pandas':pd.__version__,'h5py':h5py.__version__},'limitations':plan['limits'],'original_analysis_refits':0,'spot_level_hypothesis_tests':0}
(X/'review/analysis_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(tests.to_string(index=False),flush=True)
