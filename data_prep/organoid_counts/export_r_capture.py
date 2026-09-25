from pathlib import Path
import pickle,json,sys
import numpy as np,pandas as pd
from stream_r_object import Node,pairs,plain
def vector(n):
 if n.t==16:return plain(n)
 if n.t in [13,14,10]:
  a=pairs(n.a);v=n.v
  if isinstance(v,dict):return np.memmap(v['file'],dtype=v['dtype'],mode='r',shape=(v['length'],))
  if 'levels' in a:
   levels=plain(a['levels']);return [levels[int(i)-1] if i>0 else None for i in v]
  if n.t in [10,13]:return np.where(v==-(2**31),np.nan,v)
  return v
 raise ValueError(('Unsupported metadata vector',n.t,plain(n)))
def unpack(out):
 out=Path(out)
 with (out/'counts_parsed.pkl').open('rb') as f:c=pickle.load(f)
 with (out/'metadata_parsed.pkl').open('rb') as f:m=pickle.load(f)
 ca=pairs(c.a);ma=pairs(m.a)
 names=plain(ma['names']);obs=pd.DataFrame({k:vector(v) for k,v in zip(names,m.v)})
 obs.index=plain(ma['row.names']);obs.index.name='cell_id'
 genes,cells=plain(ca['Dimnames'])
 assert cells==list(obs.index)
 assert not obs.index.duplicated().any()
 assert len(genes)==len(set(genes))
 shape=tuple(int(x) for x in ca['Dim'].v);assert shape==(len(genes),len(cells))
 arrays={k:vector(ca[k]) for k in ['i','p','x']}
 assert len(arrays['p'])==len(cells)+1 and arrays['p'][-1]==len(arrays['i'])==len(arrays['x'])
 return arrays,obs,genes,shape
def main():
 out=Path(sys.argv[1]);arrays,obs,genes,shape=unpack(out)
 obs.to_csv(out/'cell_metadata.tsv.gz',sep='\t',compression='gzip')
 pd.DataFrame({'gene':genes}).to_csv(out/'genes.tsv',sep='\t',index=False)
 np.save(out/'column_pointers.npy',arrays['p'].astype(np.int64))
 print('SHAPE',shape,'NNZ',len(arrays['x']));print('COLUMNS',list(obs))
 for k in ['assignment','perturbation','pool','celltype','final_annotations','seurat_annotations','final.id','donor_perturbation']:
  if k in obs:print(k,obs[k].value_counts(dropna=False).to_string())
 (out/'counts_metadata_summary.json').write_text(json.dumps(dict(shape=shape,nnz=len(arrays['x']),columns=list(obs)),indent=2),encoding='utf-8')
if __name__=='__main__':main()
