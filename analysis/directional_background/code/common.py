from pathlib import Path
import numpy as np, pandas as pd
import json, hashlib, datetime, sys, os
ROOT=Path(os.environ['ICI_ANALYSIS_RUN_ROOT'])
BASE=Path(os.environ['ICI_ALL_LIGAND_RUN'])
OLD=Path(os.environ['ICI_ORGANOID_SOURCE_ROOT'])
EXT=OLD/'organoid_donor_extension'; RAW=OLD/'new_data/organoid_extracted'
TASK=Path(os.environ['ICI_SENSITIVITY_TASK_FILE'])
UP='ICI_COLITIS_EPITHELIAL_UP_STRICT';DEP='ICI_COLITIS_EPITHELIAL_LOSS_STRICT'
IFN='CURATED_IFN_VISIBILITY';RES='ICI_COLITIS_EPITHELIAL_UP_NON_IFN'
IBD1='IBD1_MEMBERSHIP_RANK';IBD2='IBD2_MEMBERSHIP_RANK'
MODS=[UP,DEP,IFN,RES,IBD1,IBD2]
SHORT=dict(zip(MODS,['Disease-up','Disease-depleted','IFN','Residual','IBD1 membership','IBD2 membership']))
FAMS=['Stem','Transit_Amplifying','Colonocyte','Goblet'];DONORS=['DNW14','DNW15','DNW9']
TGF=['BMP2','BMP3','BMP4','BMP5','BMP6','BMP8B','GDF10','GDF11','GDF15','GDF5','GDF6','INHBA','TGFB1','TGFB2','TGFB3']
SELECTED=TGF+['INS','TNFSF12','IFNG','TNFA','No_Cytomix']
def read(p):return pd.read_csv(p,sep='\t',float_precision='round_trip')
def save(d,n):
 p=ROOT/n;p.parent.mkdir(exist_ok=True,parents=True);d.to_csv(p,sep='\t',index=False,na_rep='NA',float_format='%.17g')
def jw(p,o):Path(p).write_text(json.dumps(o,indent=2,ensure_ascii=False,default=str)+'\n',encoding='utf8')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def definitions():
 genes=read(BASE/'inputs/pseudobulk_genes.tsv').gene.tolist();gm={g:i for i,g in enumerate(genes)}
 d=pd.read_csv(BASE/'inputs/frozen_module_definitions.csv');b=read(BASE/'inputs/IBD_membership_definitions.tsv')
 allsets={m:sorted(set(t.gene)) for m,t in pd.concat([d[['module','gene']],b[['module','gene']]]).groupby('module') if m in MODS}
 ix={m:np.array([gm[g] for g in allsets[m] if g in gm],dtype=np.int32) for m in MODS}
 assert len(genes)==len(set(genes)) and all(len(x)>0 for x in ix.values())
 return genes,allsets,ix
def size_factors(c):
 c=np.asarray(c,dtype=float);x=c[(c.sum(axis=1)>=10)&(c>0).all(axis=1)]
 if len(x)==0:raise ValueError('No all-positive genes for original ratio size factors')
 l=np.log(x);return np.exp(np.median(l-l.mean(axis=1,keepdims=True),axis=0))
def pattern(u,d):
 u=np.asarray(u);d=np.asarray(d)
 if len(u)!=3 or len(d)!=3 or not np.isfinite(u).all() or not np.isfinite(d).all():return 'NE'
 a=(u<0).all();b=(d>0).all()
 return 'BOTH' if a and b else 'UP_ONLY' if a else 'DEPLETED_ONLY' if b else 'OTHER_OBSERVED'
def stream_seed(family,ix):
 tag=family+'|'+','.join(map(str,sorted(ix)))
 return int.from_bytes(hashlib.sha256(tag.encode()).digest()[:8],'little')^20260826
def bg_library(expr,det,target,family,n=1000):
 # Ties share bins; sparse clusters merge before any random draw.
 eb=np.unique(np.quantile(expr,np.linspace(0,1,11)));db=np.unique(np.quantile(det,np.linspace(0,1,11)))
 e=np.searchsorted(eb[1:-1],expr,side='right');d=np.searchsorted(db[1:-1],det,side='right')
 coords=np.c_[e,d];keys=sorted(set(map(tuple,coords)));clusters={k:np.flatnonzero((coords==k).all(axis=1)) for k in keys}
 targetset=set(map(int,target));merges=[]
 while True:
  shortage=[k for k,v in clusters.items() if sum(int(g) not in targetset for g in v)<sum(int(g) in targetset for g in v)]
  if not shortage:break
  if len(clusters)==1:raise ValueError('Insufficient nonmember background genes')
  k=min(shortage);center=coords[clusters[k]].mean(axis=0)
  k2=min((j for j in clusters if j!=k),key=lambda j:(float(((coords[clusters[j]].mean(axis=0)-center)**2).sum()),j))
  merged=np.sort(np.r_[clusters.pop(k),clusters.pop(k2)]);new=min(k,k2);clusters[new]=merged
  merges.append({'from':k,'to':k2,'new_key':new})
 rng=np.random.default_rng(stream_seed(family,target));lib=np.empty((n,len(target)),dtype=np.int32);pos=0;strata=[]
 for key,v in sorted(clusters.items()):
  tar=[g for g in v if int(g) in targetset];cand=np.array([g for g in v if int(g) not in targetset],dtype=np.int32)
  if not tar:continue
  for b in range(n):lib[b,pos:pos+len(tar)]=rng.choice(cand,len(tar),replace=False)
  strata.append(dict(cluster=str(key),n_target=len(tar),n_candidates=len(cand),constituent_bins=';'.join(map(str,sorted(set(map(tuple,coords[v])))))))
  pos+=len(tar)
 assert pos==len(target) and all(len(np.unique(row))==len(target) for row in lib) and not np.isin(lib,target).any()
 return lib,strata,merges,dict(expression_breakpoints=eb.tolist(),detection_breakpoints=db.tolist())
def bg_stats(gene_delta,lib,raw):
 vals=gene_delta[lib].mean(axis=1);q=np.quantile(vals,[.025,.25,.5,.75,.975])
 return dict(B=q[2],matched_background_relative_delta=raw-q[2],random_q025=q[0],random_q25=q[1],random_q50=q[2],random_q75=q[3],random_q975=q[4],raw_reference_percentile=float(((vals<raw).sum()+.5*(vals==raw).sum())/len(vals)),random_negative_fraction=float((vals<0).mean()),random_positive_fraction=float((vals>0).mean()),random_zero_fraction=float((vals==0).mean()))
