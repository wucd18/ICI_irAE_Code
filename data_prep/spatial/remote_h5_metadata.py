"""Read public AnnData metadata via verified byte ranges, with a bounded cache."""
from pathlib import Path
import io,requests,h5py,json,hashlib,sys,numpy as np
from datetime import datetime,timezone
O=Path(__file__).resolve().parents[1]
class Remote(io.RawIOBase):
 def __init__(self,url,size,folder):
  self.url=url;self.size=size;self.pos=0;self.folder=folder;self.folder.mkdir(parents=True,exist_ok=True);self.records=[];self.session=requests.Session()
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,offset,whence=0):
  self.pos=offset if whence==0 else self.pos+offset if whence==1 else self.size+offset
  return self.pos
 def readinto(self,b):
  data=self.read(len(b));b[:len(data)]=data;return len(data)
 def read(self,n=-1):
  n=self.size-self.pos if n<0 else min(n,self.size-self.pos)
  if n<=0:return b''
  if n>150*1024**2:raise ValueError('Metadata read exceeds 150 MiB; refusing full object access')
  out=[];remain=n;block=2**20
  while remain:
   start=(self.pos//block)*block;end=min(self.size-1,start+block-1);p=self.folder/f'{start}-{end}.bin'
   if p.exists():data=p.read_bytes()
   else:
    q=self.session.get(self.url,headers={'Range':f'bytes={start}-{end}'},timeout=(20,120));q.raise_for_status()
    assert q.status_code==206 and q.headers.get('Content-Range')==f'bytes {start}-{end}/{self.size}',(q.status_code,q.headers)
    data=q.content;assert len(data)==end-start+1;p.write_bytes(data)
   self.records.append({'start':start,'end':end,'sha256':hashlib.sha256(data).hexdigest()})
   k=min(remain,len(data)-(self.pos-start));out.append(data[self.pos-start:self.pos-start+k]);self.pos+=k;remain-=k
  return b''.join(out)
def array(x):
 a=x[()]
 if a.dtype.kind in 'OS':a=np.array([v.decode() if isinstance(v,bytes) else str(v) for v in a])
 return a
def col(x):
 if isinstance(x,h5py.Dataset):return array(x)
 assert 'categories' in x and 'codes' in x,x.name
 cat=array(x['categories']);codes=array(x['codes']);return np.array([cat[k] if k>=0 else '' for k in codes])
def main():
 f=next(f for f in json.loads((O/'new_data/figshare27327813.json').read_text())['files'] if f['name']=='25_11_12_Xenium_Dataset2_5K_Annotated.h5ad')
 import pandas as pd
 folder=O/'new_data/xenium5k_metadata';folder.mkdir(exist_ok=True)
 rem=Remote(f['download_url'],f['size'],folder/'byte_ranges')
 with h5py.File(rem,'r') as h:
  print('H5 roots',list(h),flush=True)
  inventory={}
  for name in ['obs','var','X','raw','layers','obsm']:
   if name in h:
    g=h[name];inventory[name]={'attributes':{k:str(v) for k,v in g.attrs.items()},'children':list(g) if isinstance(g,h5py.Group) else [],'shape':getattr(g,'shape',None)}
  print(inventory,flush=True)
  for name in ['obs','var']:
   g=h[name];d=pd.DataFrame({k:col(g[k]) for k in g})
   d.to_csv(folder/f'{name}.tsv.gz',sep='\t',index=False)
   print(name,d.shape,[(c,d[c].nunique()) for c in d],flush=True)
  (folder/'inventory.json').write_text(json.dumps(inventory,indent=2))
 rec={'source':f,'utc':datetime.now(timezone.utc).isoformat(),'ranges':{str(v['start']):v for v in rem.records},'scope':'obs and var metadata only; no expression data or results tested'}
 (folder/'source_manifest.json').write_text(json.dumps(rec,indent=2))
if __name__=='__main__':main()
