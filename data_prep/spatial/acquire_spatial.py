"""Read GEO ZIP members over verified HTTP byte ranges; no raw reads or models."""
from pathlib import Path
import urllib.request,io,zipfile,json,time,hashlib,sys,concurrent.futures
from datetime import datetime,timezone
O=Path(__file__).resolve().parents[1]
D=O/'new_data/spatial';D.mkdir(parents=True,exist_ok=True)
class HttpResponse:
 def __init__(self,url,headers,timeout):
  with urllib.request.urlopen(urllib.request.Request(url,headers=headers),timeout=max(timeout)) as r:
   self.status_code=r.status;self.headers=dict(r.headers);self.content=r.read()
 def raise_for_status(self):
  if self.status_code>=400:raise ValueError(self.status_code)
class Session:
 def get(self,url,headers,timeout):return HttpResponse(url,headers,timeout)
class Remote(io.RawIOBase):
 def __init__(self,url):
  self.url=url;self.pos=0;self.cache={};self.session=Session()
  r=self.session.get(url,headers={'Range':'bytes=-65557'},timeout=(20,180));r.raise_for_status()
  if r.status_code!=206 or 'Content-Range' not in r.headers:raise ValueError('Range unsupported; refusing whole ZIP '+url)
  self.size=int(r.headers['Content-Range'].split('/')[-1]);self.cache[self.size-len(r.content)]=r.content
 def seekable(self):return True
 def readable(self):return True
 def tell(self):return self.pos
 def seek(self,off,whence=0):
  self.pos=off if whence==0 else self.pos+off if whence==1 else self.size+off
  return self.pos
 def read(self,n=-1):
  n=self.size-self.pos if n<0 else min(n,self.size-self.pos)
  if n<=0:return b''
  for off,data in self.cache.items():
   if off<=self.pos and off+len(data)>=self.pos+n:
    b=data[self.pos-off:self.pos-off+n];self.pos+=n;return b
  start=self.pos;end=min(self.size-1,start+max(n,65536)-1)
  for attempt in range(3):
   try:
    r=self.session.get(self.url,headers={'Range':f'bytes={start}-{end}'},timeout=(20,240));r.raise_for_status()
    assert r.status_code==206 and r.headers['Content-Range'].startswith(f'bytes {start}-{end}/')
    data=r.content;assert len(data)==end-start+1
    if len(data)<2**20:self.cache[start]=data
    self.pos+=n;return data[:n]
   except Exception:
    if attempt==2:raise
    time.sleep(2)
def one(d):
 start=datetime.now(timezone.utc).isoformat();url=d['files'][0]
 folder=D/d['gsm'];folder.mkdir(exist_ok=True)
 rec={'gsm':d['gsm'],'source_url':url,'start_utc':start,'members':[]}
 try:
  remote=Remote(url)
  with zipfile.ZipFile(remote) as z:
   inventory=[{'name':f.filename,'bytes':f.file_size,'compressed_bytes':f.compress_size,'crc32':f'{f.CRC:08x}'} for f in z.infolist()]
   rec['archive_bytes']=remote.size;rec['inventory']=inventory
   (folder/'archive_inventory.json').write_text(json.dumps(rec,indent=2),encoding='utf-8')
   if '--inventory' not in sys.argv:
    names=[f.filename for f in z.infolist() if not f.is_dir() and any(f.filename.endswith(x) for x in ['filtered_feature_bc_matrix.h5','tissue_positions_list.csv','tissue_positions.csv','scalefactors_json.json','tissue_lowres_image.png']) and not f.filename.startswith('__MACOSX/')]
    if not any(x.endswith('.h5') for x in names):
     names += [f.filename for f in z.infolist() if '/filtered_feature_bc_matrix/' in f.filename and not f.is_dir() and not f.filename.startswith('__MACOSX/')]
    assert any('matrix' in x for x in names),names
    for name in names:
     rel=Path(name)
     assert not rel.is_absolute() and '..' not in rel.parts
     dest=folder/rel;dest.parent.mkdir(parents=True,exist_ok=True)
     data=z.read(name) # zipfile verifies CRC, including all selected compressed data.
     if dest.exists() and dest.read_bytes()!=data:raise ValueError('Existing immutable download differs '+str(dest))
     if not dest.exists():dest.write_bytes(data)
     rec['members'].append({'archive_member':name,'local_path':str(dest),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'crc_verified':True})
  rec['status']='PASS';rec['exit_code']=0
 except Exception as e:
  rec['status']='FAIL';rec['error']=repr(e);rec['exit_code']=1
 rec['end_utc']=datetime.now(timezone.utc).isoformat()
 (folder/('inventory_result.json' if '--inventory' in sys.argv else 'download_result.json')).write_text(json.dumps(rec,indent=2),encoding='utf-8')
 print(d['gsm'],rec['status'],rec.get('error',''),flush=True)
 return rec
samples=json.loads((O/'review/spatial_inventory.json').read_text())
if '--one' in sys.argv:samples=[s for s in samples if s['gsm']==sys.argv[sys.argv.index('--one')+1]]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(one,samples))
(O/'review'/('spatial_zip_inventory.json' if '--inventory' in sys.argv else 'spatial_download_results.json')).write_text(json.dumps(results,indent=2),encoding='utf-8')
sys.exit(0 if all(x['status']=='PASS' for x in results) else 1)
