from common import *
from scipy.stats import rankdata
genes,full,sets=definitions();ng=len(genes)
gm=read(ROOT/'cache/pool_group_metadata.tsv');gi=gm.set_index(['donor','pool','celltype','arm'])['index'].to_dict();gc=np.load(ROOT/'cache/pool_counts.npy',mmap_mode='r')
libidx=read(ROOT/'results/matched_library_index.tsv');libs={(r.celltype,r.module):np.load(ROOT/r.path)['indices'] for r in libidx.itertuples() if r.module in MODS and r.status=='ESTIMATED'}
po=read(BASE/'inputs/all_ligand_pool_eligibility.tsv');po=po[po.eligible].copy();reuse=[]
for (d,p,f),b in po.groupby(['donor','pool','celltype']):
 reuse.append(dict(donor=d,pool=p,celltype=f,control_id='|'.join([d,p,f,'Control']),n_ligands=b.loc[b.ligand!='No_Cytomix','ligand'].nunique(),ligands=';'.join(sorted(b.loc[b.ligand!='No_Cytomix','ligand'].unique())),has_No_Cytomix='No_Cytomix' in b.ligand.values,control_cells=b.control_cells.iloc[0]))
save(pd.DataFrame(reuse),'results/control_reuse_summary.tsv')
# One rank vector for each distinct donor/pool/family/arm; duplicate controls reuse it.
rk=np.lib.format.open_memmap(ROOT/'cache/pool_ranks.npy',mode='w+',dtype=float,shape=gc.shape)
for i in range(len(gm)):rk[i,:]=rankdata(gc[i],method='average')/ng
rows=[]
for i,r in enumerate(po.itertuples()):
 ci=gi[(r.donor,r.pool,r.celltype,'Control')];ti=gi[(r.donor,r.pool,r.celltype,r.ligand)];delta=rk[ti]-rk[ci]
 for m,ix in sets.items():
  raw=delta[ix].mean();st=bg_stats(delta,libs[(r.celltype,m)],raw)
  rows.append(dict(ligand=r.ligand,celltype=r.celltype,donor=r.donor,pool=r.pool,module=m,control_id='|'.join([r.donor,r.pool,r.celltype,'Control']),treatment_cells=r.treatment_cells,control_cells=r.control_cells,treatment_score=rk[ti,ix].mean(),control_score=rk[ci,ix].mean(),delta_raw=raw,**st))
 if i%300==0:print('POOL_PAIRS',i,len(po),flush=True)
pool=pd.DataFrame(rows);save(pool,'results/pool_background_summary.tsv')
relative=[]
for keys,b in pool.groupby(['donor','celltype','pool','module']):
 for r in b.itertuples():
  ref=b[(b.ligand!=r.ligand)&(b.ligand!='No_Cytomix')];n=ref.ligand.nunique();ok=n>=3
  relative.append(dict(ligand=r.ligand,donor=keys[0],celltype=keys[1],pool=keys[2],module=keys[3],delta_raw=r.delta_raw,other_ligand_median=ref.delta_raw.median() if ok else np.nan,delta_pool_relative=r.delta_raw-ref.delta_raw.median() if ok else np.nan,n_reference_ligands=n,reference_ligands=';'.join(sorted(ref.ligand)),status='ESTIMATED' if ok else 'NE_FEWER_THAN_3_OTHER_LIGANDS'))
save(pd.DataFrame(relative),'results/pool_relative_effects.tsv')
# Independent pool contrasts and leave-one-original-pool-out, never re-fit DESeq2.
inv=read(BASE/'inputs/comparison_inventory.tsv');lp=[];lpe=[];poolclasses=[]
for r in inv.itertuples():
 b=po[(po.ligand==r.ligand)&(po.celltype==r.celltype)];ps=sorted(b.pool.unique())
 for p in ps:
  sub=pool[(pool.ligand==r.ligand)&(pool.celltype==r.celltype)&(pool.pool==p)]
  for method,col in [('raw','delta_raw'),('bg_relative','matched_background_relative_delta')]:
   w=sub.pivot(index='donor',columns='module',values=col).reindex(DONORS)
   poolclasses.append(dict(ligand=r.ligand,celltype=r.celltype,pool=p,method=method,n_matched_donors=sub.donor.nunique(),classification=pattern(w[UP],w[DEP]),reason='ESTIMATED' if sub.donor.nunique()==3 else 'NE_FEWER_THAN_3_MATCHED_DONORS'))
 for drop in ps or ['NO_ELIGIBLE_POOL']:
  left=b[b.pool!=drop];complete=left.donor.nunique()==len(DONORS)
  base=dict(ligand=r.ligand,celltype=r.celltype,dropped_pool=drop,original_eligible=r.eligibility=='ELIGIBLE',remaining_donors=';'.join(sorted(left.donor.unique())),remaining_pools=';'.join(sorted(left.pool.unique())))
  if not complete:
   lp.append(dict(**base,status='NE_FEWER_THAN_3_MATCHED_DONORS_AFTER_DELETION',raw_class='NE',bg_relative_class='NE',nonrank_class='NE'));continue
  c=[];mapping=[]
  for d in DONORS:
   poollist=sorted(left[left.donor==d].pool.unique())
   for arm in ['Control',r.ligand]:c.append(sum(gc[gi[(d,p,r.celltype,arm)]] for p in poollist));mapping.append((d,arm))
  c=np.column_stack(c);sf=size_factors(c);rr=rankdata(c,axis=0,method='average')/ng;ll=np.log2(c/sf+1);vals={m:{'raw':[],'bg_relative':[],'nonrank':[]} for m in [UP,DEP]}
  for j,d in enumerate(DONORS):
   v=rr[:,2*j+1]-rr[:,2*j];l=ll[:,2*j+1]-ll[:,2*j]
   for m,ix in sets.items():
    raw=v[ix].mean();st=bg_stats(v,libs[(r.celltype,m)],raw);non=l[ix].mean()
    lpe.append(dict(**base,donor=d,module=m,raw=raw,bg_relative=st['matched_background_relative_delta'],background_B=st['B'],nonrank=non))
    if m in vals:
     vals[m]['raw'].append(raw);vals[m]['bg_relative'].append(st['matched_background_relative_delta']);vals[m]['nonrank'].append(non)
  lp.append(dict(**base,status='ESTIMATED',**{me+'_class':pattern(vals[UP][me],vals[DEP][me]) for me in ['raw','bg_relative','nonrank']}))
save(pd.DataFrame(lp),'results/leave_one_pool_out_classification.tsv');save(pd.DataFrame(lpe),'results/leave_one_pool_out_donor_effects.tsv');save(pd.DataFrame(poolclasses),'results/independent_pool_classification.tsv')
print('POOL_COMPLETE',len(rows),'programme effects',len(lp),'deletion scenarios',flush=True)

