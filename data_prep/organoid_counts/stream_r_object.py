"""Bounded-memory XDR R serialization reader for released counts and metadata.

This reads data only; it does not evaluate any R expressions or methods. Type
layout follows R 4.4 serialize.c. Unsupported types fail explicitly. Large
numeric payloads outside the requested raw-count slots are discarded in chunks.
"""
from pathlib import Path
import struct,gzip,sys,json,pickle,hashlib,io,time
import numpy as np
from dataclasses import dataclass
sys.setrecursionlimit(10000)
@dataclass
class Node:
 t:int
 v:object=None
 a:object=None
def pairs(n):
 return {k.v:v for k,v in n.v} if n is not None else {}
def plain(n):
 if n is None:return None
 if n.t in [1,9]:return n.v
 if n.t in [13,14,10,24]:return n.v
 if n.t in [16,19,20]:return [plain(x) for x in n.v]
 if n.t in [2,6,17]:return {plain(k):plain(v) for k,v in n.v}
 if n.t==25:return {k:plain(v) for k,v in pairs(n.a).items()}
 return {'type':n.t,'value':str(n.v)}
class Captured(Exception):pass
class PartsInput(io.RawIOBase):
 def __init__(self,directory):
  self.directory=Path(directory);self.part=0;self.current=None;self.offset=0
  self.total=json.loads((self.directory/'download_progress.json').read_text(encoding='utf-8'))['total']
 def readable(self):return True
 def read(self,n=-1):
  if self.offset>=self.total:return b''
  if n<0:n=self.total-self.offset
  chunks=[];remaining=n
  while remaining and self.offset<self.total:
   if self.current is None:
    p=self.directory/f'part_{self.part:03d}.bin';start=time.time()
    while not p.exists():
     if time.time()-start>1200:raise TimeoutError(str(p))
     time.sleep(2)
    self.current=p.open('rb')
   b=self.current.read(remaining)
   if not b:self.current.close();self.current=None;self.part+=1;continue
   chunks.append(b);remaining-=len(b);self.offset+=len(b)
  return b''.join(chunks)
class Reader:
 def __init__(self,f,out,stop_after_metadata=True):
  self.f=f;self.out=Path(out);self.out.mkdir(exist_ok=True,parents=True);self.refs=[None];self.pos=0;self.nodes=0;self.counts=None;self.meta=None;self.stop_after_metadata=stop_after_metadata
  self.log=(self.out/'parse_events.jsonl').open('w',encoding='utf-8')
 def read(self,n):
  b=self.f.read(n)
  if len(b)!=n:raise EOFError((self.pos,n,len(b)))
  self.pos+=n;return b
 def integer(self):return struct.unpack('>i',self.read(4))[0]
 def length(self):
  n=self.integer()
  if n==-1:n=(self.integer()<<32)+(self.integer()&0xffffffff)
  assert n>=0,n
  return n
 def header(self):
  prefix=self.read(5);assert prefix in [b'RDX3\n',b'RDX2\n'],prefix
  assert self.read(2)==b'X\n'
  ver=self.integer();self.writer=self.integer();self.minimum=self.integer();assert ver in [2,3]
  if ver==3:self.encoding=self.read(self.integer()).decode()
 def item(self,path='',flags=None):
  off=self.pos;f=self.integer() if flags is None else flags;t=f&255;attr=bool(f&512);tag=bool(f&1024);self.nodes+=1
  if t in [254,253,252,251,250,242,241]:return None
  if t==255:
   ref=f>>8
   if ref==0:ref=self.integer()
   assert 0<ref<len(self.refs),(ref,len(self.refs),off)
   return self.refs[ref]
  if t==1:
   n=Node(t,self.item(path).v);self.refs.append(n);return n
  if t==238:
   info=self.item(path+'/altinfo');state=self.item(path+'/altstate');a=self.item(path+'/altattr')
   data=plain(info);name=list(data)[0] if isinstance(data,dict) else str(data)
   if name in ['compact_intseq','compact_realseq']:
    vals=plain(state);length,start,step=map(float,vals);n=Node(13 if name=='compact_intseq' else 14,np.arange(int(length))*step+start,a)
   else:n=Node(t,(info,state),a)
   return n
  if t in [2,6,3,5,17]:
   vals=[];a=None
   while True:
    a=self.item(path+'/attr') if attr else a
    k=self.item(path+'/tag') if tag else Node(1,'_untagged')
    sub=path+'/'+str(k.v)
    v=self.item(sub);vals.append((k,v))
    nf=self.integer();nt=nf&255
    if nt!=t:
     tail=self.item(path+'/tail',nf);assert tail is None,(nt,path,self.pos);break
    attr=bool(nf&512);tag=bool(nf&1024)
   return Node(t,vals,a)
  if t in [249,248,247]:
   assert self.integer()==0
   n=self.integer();v=[self.item(path) for _ in range(n)];r=Node(t,v);self.refs.append(r);return r
  if t==4:
   locked=self.integer();n=Node(t);self.refs.append(n);n.v=[self.item(path+'/env') for _ in range(4)];return n
  n=Node(t)
  if t in [7,8,9]:
   size=self.integer();n.v=None if size==-1 else self.read(size).decode('utf-8',errors='strict')
  elif t in [10,13,14,15,24]:
   size=self.length();dtype={10:'>i4',13:'>i4',14:'>f8',15:'>c16',24:'u1'}[t];width=np.dtype(dtype).itemsize
   self.log.write(json.dumps(dict(path=path,type=t,length=size,offset=off))+'\n');self.log.flush()
   keep='/assays/0/counts/' in path
   if size>500000:
    dest=self.out/(path.replace('/','_').replace(':','_')+'.bin') if keep else None
    sink=dest.open('xb') if dest else None;h=hashlib.sha256();remaining=size*width
    while remaining:
     b=self.read(min(8*1024*1024,remaining));remaining-=len(b)
     if sink:sink.write(b);h.update(b)
    if sink:sink.close()
    n.v=dict(file=str(dest) if dest else None,length=size,dtype=dtype,sha256=h.hexdigest() if dest else None)
   else:n.v=np.frombuffer(self.read(size*width),dtype=dtype).copy()
  elif t in [16,19,20]:
   size=self.length();n.v=[self.item(path+'/'+str(i)) for i in range(size)]
  elif t==22:
   self.refs.append(n);n.v=[self.item(path+'/pointer') for _ in range(2)]
  elif t==23:self.refs.append(n)
  elif t==25:pass
  else:raise ValueError(('Unsupported R serialized type',t,off,path))
  if attr:n.a=self.item(path)
  if path.endswith('/assays/0/counts') and t==25:
   self.counts=n
   with (self.out/'counts_parsed.pkl').open('xb') as ff:pickle.dump(n,ff)
   print('COUNTS_CAPTURED',self.pos,flush=True)
  if path.endswith('/meta.data') and t==19:
   self.meta=n
   with (self.out/'metadata_parsed.pkl').open('xb') as ff:pickle.dump(n,ff)
   print('METADATA_CAPTURED',self.pos,flush=True)
   if self.stop_after_metadata:raise Captured()
  return n
def main():
 import argparse
 p=argparse.ArgumentParser();p.add_argument('file');p.add_argument('output');a=p.parse_args()
 source=Path(a.file)
 opener=gzip.GzipFile(fileobj=PartsInput(source)) if source.is_dir() else gzip.open(source,'rb')
 with opener as f:
  r=Reader(f,a.output);r.header()
  try:r.item()
  except Captured:pass
  assert r.counts is not None and r.meta is not None
  end=r.pos
  while b:=f.read(8*1024*1024):end+=len(b)
  info=dict(status='PASS',parsed_to_metadata_bytes=r.pos,total_uncompressed_bytes=end,nodes=r.nodes,gzip_crc='PASS',writer=r.writer)
  (Path(a.output)/'parse_summary.json').write_text(json.dumps(info,indent=2),encoding='utf-8');print(info)
if __name__=='__main__':main()
