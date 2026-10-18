from common import *
import pickle,shutil
from stream_r_object import Node,pairs

genes,fullsets,sets=definitions();ng=len(genes)
for p in [TASK,BASE/'inputs/frozen_module_definitions.csv',BASE/'inputs/IBD_membership_definitions.tsv',BASE/'inputs/comparison_inventory.tsv',BASE/'inputs/all_ligand_pool_eligibility.tsv',BASE/'inputs/all_ligand_pseudobulk_metadata.tsv']:
 shutil.copy2(p,ROOT/'inputs'/p.name)
pool=read(BASE/'inputs/all_ligand_pool_eligibility.tsv');keep=pool[pool.eligible].copy()
needed=set()
for r in keep.itertuples():
 needed.add((r.donor,r.pool,r.celltype,r.ligand));needed.add((r.donor,r.pool,r.celltype,'Control'))
usecols=['cell_id','assignment','pool','celltype','final_annotations','perturbation']+[f'X{i}_DEG' for i in range(1,21)]
obs=pd.read_csv(RAW/'cell_metadata.tsv.gz',sep='\t',usecols=usecols)
assert not obs.cell_id.duplicated().any() and obs[usecols].notna().all().all()
obs['arm']=obs.perturbation.where(~obs.perturbation.isin(['Sato_1','Sato_2','Sato_3']),'Control')
cnt=obs.groupby(['assignment','pool','celltype','final_annotations','arm'],sort=True).size().reset_index(name='n_cells')
save(cnt,'results/recovered_substate_counts.tsv')
state_mapping=obs[['celltype','final_annotations']].drop_duplicates().sort_values(['celltype','final_annotations']);save(state_mapping,'inputs/released_annotation_mapping.tsv')
with (RAW/'counts_parsed.pkl').open('rb') as f:parsed=pickle.load(f)
attrs=pairs(parsed.a)
def mem(n):
 v=attrs[n].v;return np.memmap(v['file'],dtype=v['dtype'],mode='r',shape=(v['length'],))
ii=mem('i');xx=mem('x');ptr=np.load(RAW/'column_pointers.npy');assert len(ptr)==len(obs)+1
keys=sorted(needed);ki={k:i for i,k in enumerate(keys)}
gc=np.lib.format.open_memmap(ROOT/'cache/pool_counts.npy',mode='w+',dtype=np.int64,shape=(len(keys),ng))
met=[];detect={};statelog=[];gp=[]
for k,positions in obs.groupby(['assignment','pool','celltype','arm'],sort=True).indices.items():
 if k not in needed:continue
 ix=ki[k];sub=obs.iloc[positions];gc[ix,:]=0;iscontrol=k[3]=='Control';comp=k[3] in SELECTED or iscontrol
 dc=np.zeros(ng,dtype=np.int32) if iscontrol else None
 ls={state:np.zeros(ng,dtype=float) for state in sub.final_annotations.unique()} if comp else {}
 for cell in positions:
  start,end=int(ptr[cell]),int(ptr[cell+1]);gi=ii[start:end];v=xx[start:end]
  assert (v>0).all() and np.isfinite(v).all() and np.equal(v,np.floor(v)).all()
  gc[ix,gi]+=v.astype(np.int64)
  if iscontrol:dc[gi]+=1
  if comp:ls[obs.final_annotations.iat[cell]][gi]+=np.log1p(10000*v/v.sum())
 met.append(dict(index=ix,donor=k[0],pool=k[1],celltype=k[2],arm=k[3],n_cells=len(positions),control_id='|'.join(k[:3])+'|Control'))
 if iscontrol:detect[ix]=dc
 if comp:
  for state,b in sub.groupby('final_annotations'):
   for m,g in sets.items():statelog.append(dict(donor=k[0],pool=k[1],celltype=k[2],arm=k[3],substate=state,module=m,n_cells=len(b),mean_cell_log_normalized=ls[state][g].sum()/len(g)/len(b)))
   for n in range(1,21):gp.append(dict(donor=k[0],pool=k[1],celltype=k[2],arm=k[3],substate=state,programme=f'GP{n}',n_cells=len(b),saved_cNMF_usage=float(b[f'X{n}_DEG'].mean())))
 if len(met)%100==0:print('POOL_GROUPS',len(met),len(keys),flush=True)
assert len(met)==len(keys)
gm=pd.DataFrame(met).sort_values('index');save(gm,'cache/pool_group_metadata.tsv');gc.flush()
save(pd.DataFrame(statelog),'cache/substate_log_expression.tsv');save(pd.DataFrame(gp),'cache/substate_saved_programmes.tsv')
# Unique matched controls from eligible three-donor comparisons only.
inv=read(BASE/'inputs/comparison_inventory.tsv');valid=inv[inv.eligibility.eq('ELIGIBLE')][['ligand','celltype']]
validpool=keep.merge(valid,on=['ligand','celltype']);brows=[];baselines={}
for fam in FAMS:
 ex=[];de=[]
 for donor in DONORS:
  ps=sorted(validpool[(validpool.celltype==fam)&(validpool.donor==donor)].pool.unique());idx=[ki[(donor,p,fam,'Control')] for p in ps]
  c=np.asarray(gc[idx,:]).sum(axis=0);nc=int(gm.set_index('index').loc[idx,'n_cells'].sum());det=sum(detect[i] for i in idx)/nc
  expr=np.log2(1+1e6*c/c.sum());ex.append(expr);de.append(det)
  brows.append(dict(celltype=fam,donor=donor,pools=';'.join(ps),n_unique_controls=len(idx),n_cells=nc,control_ids=';'.join('|'.join((donor,p,fam,'Control')) for p in ps)))
 baselines[fam+'_expression']=np.mean(ex,axis=0);baselines[fam+'_detection']=np.mean(de,axis=0)
np.savez_compressed(ROOT/'cache/matching_baselines.npz',**baselines);save(pd.DataFrame(brows),'results/matching_baseline_controls.tsv')
np.savez_compressed(ROOT/'cache/control_detection_counts.npz',**{str(k):v for k,v in detect.items()})
print('PREPARATION_COMPLETE',len(keys),'groups',len(statelog),'substate scores',flush=True)
