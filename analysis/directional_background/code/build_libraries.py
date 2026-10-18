from common import *
genes,full,sets=definitions();baseline=np.load(ROOT/'cache/matching_baselines.npz')
allsets=dict(sets);pairs=[]
for a in [UP,DEP]:
 for b in [IFN,IBD1,IBD2]:
  an=a+'__without__'+b;bn=b+'__without__'+a
  allsets[an]=np.setdiff1d(sets[a],sets[b]);allsets[bn]=np.setdiff1d(sets[b],sets[a])
  pairs.append(dict(left=a,right=b,left_exclusive=an,right_exclusive=bn,n_left_remaining=len(allsets[an]),n_right_remaining=len(allsets[bn])))
save(pd.DataFrame(pairs),'inputs/overlap_sensitivity_pairs.tsv')
save(pd.DataFrame({'index':np.arange(len(genes)),'gene':genes}),'inputs/gene_universe.tsv')
index=[];balance=[];strata=[];mergelog=[];breakpoints={}
(ROOT/'results/matched_set_membership').mkdir(exist_ok=True)
for fam in FAMS:
 expr=baseline[fam+'_expression'];det=baseline[fam+'_detection'];cache={}
 for name,ix in allsets.items():
  key=hashlib.sha256(np.sort(ix).astype('<i4').tobytes()).hexdigest()[:16]
  status='ESTIMATED' if len(ix)>=10 else 'NE_FEWER_THAN_10_REMAINING_MEMBERS'
  f=ROOT/'results/matched_set_membership'/f'{fam}_{key}.npz'
  index.append(dict(celltype=fam,module=name,library_id=f'{fam}_{key}',path=str(f.relative_to(ROOT)),n_members=len(ix),n_sets=1000 if status=='ESTIMATED' else 0,status=status))
  if status!='ESTIMATED':continue
  if key not in cache:
   lib,st,merges,bp=bg_library(expr,det,ix,fam);np.savez_compressed(f,indices=lib,target_indices=ix)
   cache[key]=lib
   for s in st:strata.append(dict(celltype=fam,library_id=f'{fam}_{key}',**s))
   for i,m in enumerate(merges):mergelog.append(dict(celltype=fam,library_id=f'{fam}_{key}',merge_order=i,**m))
   breakpoints[fam]=bp
  lib=cache[key]
  for dim,v in [('log2_control_CPM_plus1',expr),('cell_detection_rate',det)]:
   means=v[lib].mean(axis=1);tv=v[ix]
   balance.append(dict(celltype=fam,module=name,library_id=f'{fam}_{key}',dimension=dim,n_target=len(ix),target_mean=tv.mean(),background_mean=means.mean(),background_mean_q025=np.quantile(means,.025),background_mean_q975=np.quantile(means,.975),target_sd=tv.std(),mean_difference=means.mean()-tv.mean(),standardized_mean_difference=(means.mean()-tv.mean())/v.std() if v.std()>0 else np.nan,target_q10=np.quantile(tv,.1),target_q90=np.quantile(tv,.9),background_member_q10=np.quantile(v[lib],.1),background_member_q90=np.quantile(v[lib],.9)))
 print('MATCH_LIBRARIES',fam,len(cache),flush=True)
save(pd.DataFrame(index),'results/matched_library_index.tsv');save(pd.DataFrame(balance),'results/matching_balance.tsv');save(pd.DataFrame(strata),'results/matching_strata.tsv');save(pd.DataFrame(mergelog,columns=['celltype','library_id','merge_order','from','to','new_key']),'results/matching_sparse_bin_merges.tsv');jw(ROOT/'results/matching_breakpoints.json',breakpoints)

