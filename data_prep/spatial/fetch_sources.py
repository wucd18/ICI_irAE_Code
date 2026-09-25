from pathlib import Path
import requests,concurrent.futures,json,hashlib,sys
from datetime import datetime,timezone
O=Path(__file__).resolve().parents[1]
jobs=[
 ('mendeley_files.json','https://data.mendeley.com/public-api/datasets/7z8yx644hb/files?folder_id=root&version=1'),
 ('GSM4797916_A1.tar.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4797nnn/GSM4797916/suppl/GSM4797916_A1.tar.gz'),
 ('GSM4797917_A2.tar.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4797nnn/GSM4797917/suppl/GSM4797917_A2.tar.gz')]
if len(sys.argv)>1:jobs=json.loads(Path(sys.argv[1]).read_text())
def one(job):
 name,url=job;rec={'url':url,'file':name,'start_utc':datetime.now(timezone.utc).isoformat()}
 try:
  r=requests.get(url,timeout=(20,240));r.raise_for_status();p=O/'new_data'/name
  if p.exists() and p.read_bytes()!=r.content:raise ValueError('Existing file differs; new filename required')
  if not p.exists():p.write_bytes(r.content)
  rec.update(status='PASS',bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest(),http_status=r.status_code)
 except Exception as e:rec.update(status='FAIL',error=repr(e))
 rec['end_utc']=datetime.now(timezone.utc).isoformat();print(name,rec['status'],rec.get('error',''),flush=True);return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:out=list(pool.map(one,jobs))
(O/'logs'/('fetch_sources_'+datetime.now(timezone.utc).strftime('%H%M%S%f')+'.json')).write_text(json.dumps(out,indent=2),encoding='utf-8')
sys.exit(0 if all(r['status']=='PASS' for r in out) else 1)
