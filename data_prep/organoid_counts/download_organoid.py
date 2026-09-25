from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import requests,hashlib,json,time,os,sys
from datetime import datetime,timezone
sys.stdout.reconfigure(encoding='utf-8')
O=Path(__file__).resolve().parents[1];D=O/'new_data/organoid_download';D.mkdir(exist_ok=True)
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313368/suppl/GSE313368_Combined_Pool_Merge_Harmonized_Donor_Only.Robj.gz'
probe=requests.get(url,headers={'Range':'bytes=0-0','Accept-Encoding':'identity'},timeout=45)
assert probe.status_code==206 and probe.headers['Content-Range'].startswith('bytes 0-0/')
total=int(probe.headers['Content-Range'].split('/')[-1]);chunk=128*1024*1024
spec=[(i,s,min(total-1,s+chunk-1)) for i,s in enumerate(range(0,total,chunk))]
def fetch(t):
 i,s,e=t;dest=D/f'part_{i:03d}.bin'
 if dest.exists() and dest.stat().st_size==e-s+1:return dict(part=i,start=s,end=e,size=e-s+1,reused=True)
 errors=[]
 for attempt in range(1,4):
  p=D/f'part_{i:03d}_attempt_{attempt}.partial'
  if p.exists():continue
  try:
   with requests.get(url,headers={'Range':f'bytes={s}-{e}','Accept-Encoding':'identity'},stream=True,timeout=(30,120)) as r:
    r.raise_for_status();assert r.status_code==206 and r.headers['Content-Range']==f'bytes {s}-{e}/{total}'
    h=hashlib.sha256()
    with p.open('xb') as f:
     for b in r.iter_content(1024*1024):f.write(b);h.update(b)
   assert p.stat().st_size==e-s+1
   p.rename(dest);return dict(part=i,start=s,end=e,size=e-s+1,sha256=h.hexdigest(),attempt=attempt,prior_errors=errors)
  except Exception as err:errors.append(repr(err))
 raise RuntimeError((i,errors))
started=time.time();rows=[]
with ThreadPoolExecutor(max_workers=4) as ex:
 for f in as_completed([ex.submit(fetch,t) for t in spec]):
  rows.append(f.result());done=sum(x['size'] for x in rows);print(f'{done}/{total} bytes; {done/max(1,time.time()-started)/1e6:.1f} MB/s',flush=True)
  (D/'download_progress.json').write_text(json.dumps(dict(url=url,total=total,completed_bytes=done,parts=rows),indent=2),encoding='utf-8')
dest=O/'new_data/GSE313368_Combined_Pool_Merge_Harmonized_Donor_Only.Robj.gz'
h=hashlib.sha256()
with dest.open('xb') as f:
 for i,s,e in spec:
  with (D/f'part_{i:03d}.bin').open('rb') as inp:
   while b:=inp.read(8*1024*1024):f.write(b);h.update(b)
assert dest.stat().st_size==total
(D/'download_complete.json').write_text(json.dumps(dict(file=str(dest),url=url,size=total,sha256=h.hexdigest(),utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.time()-started),indent=2),encoding='utf-8')
print('DOWNLOAD_PASS',total,h.hexdigest(),flush=True)
