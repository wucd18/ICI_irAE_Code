from common import *
from scipy.stats import rankdata
from scipy.io import mmread
import gzip
genes,full,sets=definitions();ng=len(genes)
inv=read(BASE/'inputs/comparison_inventory.tsv');eligible=inv[inv.eligibility.eq('ELIGIBLE')]
oldmeta=read(EXT/'source_data/pseudobulk_metadata.tsv');oldc=mmread(EXT/'source_data/pseudobulk_counts.mtx').tocsc()
jobs=read(BASE/'inputs/new_model_jobs.tsv').set_index(['ligand','celltype'])
oldi=read(BASE/'inputs/targeted_model_inventory.tsv');newi=read(BASE/'results/new_model_inventory.tsv');sfref=pd.concat([oldi,newi]).set_index(['ligand','celltype'])
libindex=read(ROOT/'results/matched_library_index.tsv');libs={};allsets={}
for r in libindex.itertuples():
 if r.status!='ESTIMATED':continue
 z=np.load(ROOT/r.path);libs[(r.celltype,r.module)]=z['indices'];allsets[r.module]=z['target_indices']
raw=read(BASE/'results/all_readout_donor_deltas_long.tsv');rawindex=raw.set_index(['ligand','celltype','donor','module'])
gm=read(ROOT/'cache/pool_group_metadata.tsv');gindex=gm.set_index(['donor','pool','celltype','arm'])['index'].to_dict();gc=np.load(ROOT/'cache/pool_counts.npy',mmap_mode='r')
n=len(eligible)*len(DONORS);rmat=np.lib.format.open_memmap(ROOT/'cache/paired_gene_rank_deltas.npy',mode='w+',dtype=float,shape=(n,ng));lmat=np.lib.format.open_memmap(ROOT/'cache/paired_gene_log_deltas.npy',mode='w+',dtype=float,shape=(n,ng))
rows=[];bg=[];non=[];exclusive=[];sfrows=[];checks=[];counter=0;header=True
memberpath=ROOT/'results/member_gene_log_changes.tsv'
with memberpath.open('w',encoding='utf8',newline='') as memberfile:
 for i,r in enumerate(eligible.itertuples()):
  key=(r.ligand,r.celltype)
  if r.action=='REUSE_TARGETED':
   b=oldmeta[(oldmeta.ligand==r.ligand)&(oldmeta.celltype==r.celltype)].copy();c=oldc[:,b.index].toarray()
  else:
   j=jobs.loc[key];b=read(j.metadata);df=read(j.counts);assert df.gene.tolist()==genes and df.columns[1:].tolist()==b.sample_id.tolist();c=df.iloc[:,1:].to_numpy(dtype=np.int64)
  assert len(b)==len(DONORS)*2 and not b.duplicated(['donor','treatment']).any()
  pooled=np.column_stack([sum(gc[gindex[(x.donor,p,x.celltype,'Control' if x.treatment=='Control' else x.ligand)]] for p in x.pools.split(';')) for x in b.itertuples()])
  assert np.array_equal(c,pooled),key
  sf=size_factors(c);norm=c/sf;logs=np.log2(norm+1);ranks=rankdata(c,axis=0,method='average')/ng
  ref=sfref.loc[key];err=max(abs(sf.min()-ref.min_size_factor),abs(sf.max()-ref.max_size_factor))
  assert err<1e-11,(key,err)
  checks.append(dict(ligand=r.ligand,celltype=r.celltype,pool_counts_equal_saved=True,saved_factor_extrema_max_abs_error=err))
  for j,x in enumerate(b.itertuples()):sfrows.append(dict(ligand=r.ligand,celltype=r.celltype,donor=x.donor,treatment=x.treatment,sample_id=x.sample_id,pools=x.pools,size_factor=sf[j],normalization='original_DESeq2_ratio_formula_only_no_model_refit'))
  for donor in DONORS:
   ci=np.flatnonzero((b.donor==donor)&(b.treatment=='Control'))[0];ti=np.flatnonzero((b.donor==donor)&(b.treatment=='Treatment'))[0]
   x=b.iloc[ti];pools=x.pools;controlids=';'.join('|'.join((donor,p,r.celltype,'Control')) for p in pools.split(';'))
   rankdelta=ranks[:,ti]-ranks[:,ci];logdelta=logs[:,ti]-logs[:,ci];rmat[counter,:]=rankdelta;lmat[counter,:]=logdelta
   base=dict(pair_index=counter,ligand=r.ligand,celltype=r.celltype,donor=donor,pools=pools,control_ids=controlids,condition_class=r.condition_class)
   rows.append(base)
   for m,ix in sets.items():
    original=rawindex.loc[(r.ligand,r.celltype,donor,m)];rv=float(original.delta);assert abs(rankdelta[ix].mean()-rv)<2e-14,(key,donor,m)
    st=bg_stats(rankdelta,libs[(r.celltype,m)],rv);bg.append(dict(**base,module=m,delta_raw=rv,treatment_score=original.treatment_score,control_score=original.control_score,**st,n_members=len(ix)))
    v=logdelta[ix];lo=(norm[ix,ci]<1)&(norm[ix,ti]<1)
    non.append(dict(**base,module=m,normalized_expression_change=v.mean(),member_median=v.mean() if len(v)==1 else np.median(v),member_q025=np.quantile(v,.025),member_q25=np.quantile(v,.25),member_q75=np.quantile(v,.75),member_q975=np.quantile(v,.975),positive_fraction=(v>0).mean(),negative_fraction=(v<0).mean(),zero_fraction=(v==0).mean(),n_members=len(ix),both_arms_normalized_count_below1_fraction=lo.mean(),treatment_size_factor=sf[ti],control_size_factor=sf[ci]))
    tab=pd.DataFrame(dict(ligand=r.ligand,celltype=r.celltype,donor=donor,module=m,gene=np.array(genes)[ix],control_normalized_count=norm[ix,ci],treatment_normalized_count=norm[ix,ti],gene_log_delta=v,both_arms_below1=lo))
    tab.to_csv(memberfile,sep='\t',index=False,header=header,float_format='%.17g');header=False
   for m,ix in allsets.items():
    if m in MODS:continue
    v=float(rankdelta[ix].mean());st=bg_stats(rankdelta,libs[(r.celltype,m)],v)
    exclusive.append(dict(**base,module=m,n_members=len(ix),raw=v,bg_relative=st['matched_background_relative_delta'],background_B=st['B'],nonrank=logdelta[ix].mean()))
   counter+=1
  if i%30==0:print('SCORED',i,'/',len(eligible),flush=True)
assert counter==n
rmat.flush();lmat.flush();save(pd.DataFrame(rows),'cache/paired_index.tsv');save(pd.DataFrame(bg),'results/matched_background_effects.tsv');save(pd.DataFrame(non),'results/nonrank_donor_effects.tsv');save(pd.DataFrame(exclusive),'results/overlap_removed_donor_effects.tsv');save(pd.DataFrame(sfrows),'results/normalization_factors.tsv');save(pd.DataFrame(checks),'review/paired_input_and_normalization_checks.tsv')
save(raw,'results/raw_donor_effects.tsv')
print('MAIN_SCORES_COMPLETE',n,'pairs',flush=True)
