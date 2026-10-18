from pathlib import Path
import json,hashlib,datetime,os,sys
import numpy as np,pandas as pd
ROOT=Path(os.environ['ICI_ANALYSIS_RUN_ROOT'])
PROJECT=Path(os.environ['ICI_PROJECT_ROOT'])
OLD=Path(os.environ['ICI_ORGANOID_SOURCE_ROOT'])
E=OLD/'organoid_donor_extension';SAVED=E/'results/20260914T062958'
RAW=OLD/'new_data/organoid_extracted'
CODEBASE=Path(os.environ['ICI_RECORDED_CODE_ROOT'])
TASK=Path(os.environ['ICI_ANALYSIS_TASK_ROOT'])
FROZEN=Path(os.environ['ICI_FROZEN_MEMBERSHIP'])
for n in ['inputs','results','figures','figure_source_data','text','review','logs','temp','inputs/model_counts','results/new_models']:(ROOT/n).mkdir(parents=True,exist_ok=True)
def jw(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=str)+'\n',encoding='utf8')
def jr(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def save(df,name):df.to_csv(ROOT/name,sep='\t',index=False,na_rep='NA',float_format='%.17g')
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def register(paths):
 out=ROOT/'review/input_manifest.tsv'
 records=pd.read_csv(out,sep='\t').to_dict('records') if out.exists() else []
 known={r['path'] for r in records}
 for p,role in paths:
  p=Path(p);assert p.is_file(),p
  if str(p) not in known:
   records.append(dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p),source='GSE313368' if 'organoid' in str(p) or 'GSE313368' in str(p) else 'approved project input',role=role));known.add(str(p))
 save(pd.DataFrame(records),'review/input_manifest.tsv')

