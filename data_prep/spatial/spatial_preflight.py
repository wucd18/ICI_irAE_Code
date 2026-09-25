import os
"""Freeze complete sample mapping and scoring scope before viewing module outcomes."""
from pathlib import Path
import json,tarfile,hashlib,gzip,re
import pandas as pd,h5py,numpy as np
from datetime import datetime,timezone
O=Path(__file__).resolve().parents[1];S=Path(os.environ['ICI_ARCHIVE_ROOT'])/'original_v1'
X=O/'spatial_extension';X.mkdir(exist_ok=True)
for d in ['source_data','figures','review']: (X/d).mkdir(exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
rows=json.loads((O/'review/publisher_spatial_samples.json').read_text())
geo=json.loads((O/'review/spatial_inventory.json').read_text());bysect={r['section']:r for r in geo}
inputs=[];meta=[]
for gsm,sec in [('GSM4797916','A1'),('GSM4797917','A2')]:
 archive=O/'new_data'/f'{gsm}_{sec}.tar.gz';dest=O/'new_data/spatial'/gsm;dest.mkdir(exist_ok=True)
 with tarfile.open(archive) as tar:
  for m in tar.getmembers():
   if not m.isfile() or not any(m.name.endswith(x) for x in ['matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz','tissue_positions_list.csv','scalefactors_json.json','tissue_lowres_image.png']):continue
   p=dest/m.name
   assert p.resolve().is_relative_to(dest.resolve())
   data=tar.extractfile(m).read();p.parent.mkdir(parents=True,exist_ok=True)
   if p.exists():assert p.read_bytes()==data
   else:p.write_bytes(data)
 bysect[sec]={'gsm':gsm,'section':sec,'characteristics':{'disease':'Healthy'},'files':[f'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4797nnn/{gsm}/suppl/{gsm}_{sec}.tar.gz']}
for row in rows:
 for sec in row['Hashing Code'].replace(' ','').split('&'):
  g=bysect[sec];folder=O/'new_data/spatial'/g['gsm']/sec
  h5=folder/'filtered_feature_bc_matrix.h5'
  if h5.exists():
   with h5py.File(h5) as f:
    ids=[x.decode() for x in f['matrix/features/id'][:]];sym=[x.decode() for x in f['matrix/features/name'][:]];shape=f['matrix/shape'][:].tolist()
   matrix=str(h5);features=barcodes='';source_processing='released Space Ranger filtered matrix'
  else:
   features=str(folder/'raw_feature_bc_matrix/features.tsv.gz');barcodes=str(folder/'raw_feature_bc_matrix/barcodes.tsv.gz');matrix=str(folder/'raw_feature_bc_matrix/matrix.mtx.gz')
   f=pd.read_csv(features,sep='\t',header=None);ids=f[0].tolist();sym=f[1].tolist();shape=[len(ids),len(pd.read_csv(barcodes,header=None))];source_processing='released raw matrix, tissue flag applied'
  assert len(set(ids))==len(ids),'Duplicate feature IDs'
  condition={'CPI':'ICI colitis','Healthy':'Healthy','UC':'UC'}[row['Disease'].strip()]
  assert {'Checkpoint inhibitor induced colitis':'ICI colitis','Healthy':'Healthy','Healthy Control':'Healthy','Ulcerative Colitis':'UC'}[g['characteristics']['disease']]==condition
  coords=str(folder/'spatial/tissue_positions_list.csv');scale=str(folder/'spatial/scalefactors_json.json');image=str(folder/'spatial/tissue_lowres_image.png')
  meta.append(dict(section=sec,gsm=g['gsm'],patient=row['Study Code'],condition=condition,age=row['Age'],sex=row['Gender'].strip(),source_series='GSE158328' if sec.startswith('A') else 'GSE189184',source_url=g['files'][0],source_excel_row=row['source_excel_row'],source_sheet='Patient Data',source_workbook='new_data/Gupta_Table_S1_publisher.xlsx',matrix=matrix,features=features,barcodes=barcodes,coordinates=coords,scalefactors=scale,image=image,processing=source_processing,n_features=shape[0],n_matrix_spots=shape[1],duplicate_symbols=len(sym)-len(set(sym))))
  for p in [matrix,features,barcodes,coords,scale,image]:
   if p:inputs.append(dict(path=p,bytes=Path(p).stat().st_size,sha256=sha(Path(p)),accession=g['gsm'],source_url=g['files'][0],role='frozen programme spatial extension input'))
  pd.DataFrame({'feature_id':ids,'symbol':sym,'feature_row':np.arange(len(ids))+1}).to_csv(X/'source_data'/f'{sec}_feature_index.tsv',sep='\t',index=False)
meta=pd.DataFrame(meta).sort_values('section',key=lambda x:x.str[0]+x.str[1:].str.zfill(2));assert not meta.section.duplicated().any()
meta.to_csv(X/'source_data/sample_manifest.tsv',sep='\t',index=False)
modules=Path(os.environ['ICI_FROZEN_MODULES']);d=pd.read_csv(modules,sep='\t');names=['CURATED_IFN_VISIBILITY','ICI_COLITIS_EPITHELIAL_UP_STRICT','ICI_COLITIS_EPITHELIAL_UP_NON_IFN','ICI_COLITIS_EPITHELIAL_LOSS_STRICT']
definition=d[d.module.isin(names)].copy();definition['original_source_file']=str(modules);definition['original_source_row']=definition.index+1
definition.to_csv(X/'source_data/frozen_module_members.tsv',sep='\t',index=False)
for p in [modules,O/'new_data/Gupta_Table_S1_publisher.xlsx',Path(os.environ['ICI_ORIGINAL_CODE_ROOT'])/'00_scripts/evaluate_patient_level_module_performance.py',Path(os.environ['ICI_ORIGINAL_CODE_ROOT'])/'00_scripts/project_config.R']:
 inputs.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),accession='V1 or source publication',source_url='original file or publisher supplementary table',role='frozen definition, scoring specification or patient mapping'))
pd.DataFrame(inputs).drop_duplicates('path').to_csv(X/'review/input_manifest_before.tsv',sep='\t',index=False)
plan={'locked_utc':datetime.now(timezone.utc).isoformat(),'status':'PRE_OUTCOME_LOCKED','modules':names,'sample_counts':meta.groupby('condition').size().to_dict(),'patient_counts':meta.groupby('condition').patient.nunique().to_dict(),'primary_comparison':'ICI colitis versus Healthy; UC shown as separate context without a new UC hypothesis test','replicate':'patient; sum counts across all tissue spots and both repeated sections before scoring','score':'same V1 normalized within-sample average gene rank; logCPM monotone equivalent to counts; ties use average rank; disease-depleted score=1-raw','gene_universe':'all common Ensembl IDs across released feature lists, before observing module values; no expression or outcome gene filtering','spot_inclusion':'released in_tissue=1 spots with positive library; no new clustering, DE, model fitting or histological annotation; source post-histology exclusions unavailable','coverage':'at least 10 module genes, original minimum; all matching members retained','test':'V1 two-sided scipy.stats.mannwhitneyu default method; BH across the four predeclared module comparisons in this new cohort; original statistics not changed; no bootstrap','AUC':'not generated; descriptive rank score and patient differences only','minimum_patients_per_group':3,'spatial_inference':'spot scores only visual; no spot-level P value, no epithelial ROI inferred from tested module','limits':['small cohort, healthy controls not ICI-exposed','sex and batch imbalance; additional healthy donor A1/A2 from original source study GSE158328','whole-tissue score cannot distinguish state from cell composition','public portal annotation download unresolved; no claim of epithelial-specific localization','source studies Oxford, UCSF discovery, MGH external are distinct; no common patient IDs; genomic identity linkage not available']}
(X/'review/analysis_plan_locked.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
print(json.dumps(plan,indent=2))
