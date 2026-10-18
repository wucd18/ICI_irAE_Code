from common import *
import shutil,itertools,pickle
from scipy.stats import rankdata
from scipy.io import mmread
from stream_r_object import Node,pairs

plan=jr(E/'review/analysis_plan_locked.json');families=plan['celltypes'];modules=plan['modules']
primary=['ICI_COLITIS_EPITHELIAL_UP_STRICT','ICI_COLITIS_EPITHELIAL_LOSS_STRICT','CURATED_IFN_VISIBILITY']
group=pd.read_csv(E/'source_data/all_released_group_cell_counts.tsv',sep='\t')
ligands=sorted(set(group.arm)-{'Control','No_Cytomix'})
assert len(ligands)==79,'Released condition universe differs from task; stop rather than invent'
cfg=dict(created_utc=now(),reuse_authorization='Explicit user request 2026-10-16; reuse reviewed data, code and outputs',
 prior_plan=plan,ligands=ligands,reference='No_Cytomix',families=families,modules_in_original_fgsea_order=modules,primary_modules=primary,
 FDR=dict(main='All released ligand x four family x three primary readouts; BH with full planned n, including NE',main_n=len(ligands)*len(families)*len(primary),residual='Separate all-ligand residual family',residual_n=len(ligands)*len(families),assay_reference='Separate No_Cytomix four-readout x family family, same as targeted',reference_n=len(families)*len(modules)),
 classification='Strict unrounded donor deltas; require 3 observed pairs; both if all three same donor records have up<0 and depleted>0. Otherwise up-only, depleted-only, other observed directions. Zero counts in neither sign.',
 benchmark='Within family across ligand mean paired raw-rank deltas, Pearson and Spearman descriptive coefficients without independence-based P; donor-stratified coefficients also reported. NoCytomix excluded.',
 example_selection='Before analysis: within each family and joint-pattern category choose lexicographically first ligand, plus lexicographically first ligand for each disease/IFN donor-direction discordance category. Never rank by efficacy or significance.',
 new_score_tests='None. No permutation/sign-flip analysis is a completion requirement.',old_models='Reuse all saved targeted estimates/scores/coverage/gene results unchanged',old_FDR='Keep separate immutable columns and archived tables',
 IBD_policy='Use only saved reproducible members plus required directions/weights; publisher aggregate scores are a separate source and not same-rule comparisons',run_all=False)
jw(ROOT/'analysis_config.json',cfg)
shutil.copy2(FROZEN,ROOT/'inputs/frozen_module_definitions.csv')
for n in ['module_GSEA.tsv','donor_module_scores.tsv','donor_paired_effects.tsv','module_coverage.tsv','model_inventory.tsv']:
 shutil.copy2(SAVED/n,ROOT/'inputs'/('targeted_'+n))
for n in ['Codex_Incremental_Analysis_and_Figures.md','Manuscript_Text_and_Figure4_Revision.md']:shutil.copy2(TASK/n,ROOT/'inputs'/n)
shutil.copy2(TASK/'csbj_20261016/Manuscript_CSBJ_revised_clean.docx',ROOT/'inputs/Manuscript_CSBJ_revised_clean.docx')
definitions=pd.read_csv(FROZEN);genes=pd.read_csv(RAW/'genes.tsv',sep='\t').gene.tolist();assert len(genes)==len(set(genes))
genemap={g:i for i,g in enumerate(genes)};sets={m:[genemap[g] for g in definitions.loc[definitions.module==m,'gene'].drop_duplicates() if g in genemap] for m in modules}
save(pd.DataFrame({'gene':genes}),'inputs/pseudobulk_genes.tsv')
oldmeta=pd.read_csv(E/'source_data/pseudobulk_metadata.tsv',sep='\t');oldcounts=mmread(E/'source_data/pseudobulk_counts.mtx').tocsc()
oldgenes=pd.read_csv(E/'source_data/pseudobulk_genes.tsv',sep='\t').gene.tolist();assert oldgenes==genes
oldmodel=pd.read_csv(SAVED/'model_inventory.tsv',sep='\t');oldkeys=set(zip(oldmodel.ligand,oldmodel.celltype))
oldscores=pd.read_csv(SAVED/'donor_module_scores.tsv',sep='\t')
donors=sorted(group.assignment.unique());pools=sorted(group.pool.unique());counts=group.set_index(['assignment','pool','celltype','arm']).n_cells.to_dict()
eligibility=[];comparisons=[];sample_specs=[];needed_groups=set()
for ligand,ct in itertools.product(ligands+['No_Cytomix'],families):
 paired=[]
 for donor in donors:
  keep=[]
  for pool in pools:
   nt=int(counts.get((donor,pool,ct,ligand),0));nc=int(counts.get((donor,pool,ct,'Control'),0))
   if nt:
    ok=nt>=20 and nc>=20
    eligibility.append(dict(ligand=ligand,celltype=ct,donor=donor,pool=pool,treatment_cells=nt,control_cells=nc,eligible=ok,reason='PASS' if ok else 'one_arm_below_existing_20_cell_minimum'))
    if ok:keep.append(pool)
  if keep:
   paired.append(donor)
   for arm in ['Control',ligand]:
    spec=dict(ligand=ligand,celltype=ct,donor=donor,treatment='Control' if arm=='Control' else 'Treatment',pools=';'.join(keep),n_cells=sum(counts[(donor,p,ct,arm)] for p in keep))
    sample_specs.append(spec)
    if not ((oldmeta.ligand==ligand)&(oldmeta.celltype==ct)&(oldmeta.donor==donor)&(oldmeta.treatment==spec['treatment'])).any():
     needed_groups.update((donor,p,ct,arm) for p in keep)
 comparisons.append(dict(ligand=ligand,condition_class='assay_reference' if ligand=='No_Cytomix' else 'ligand',celltype=ct,n_donors=len(paired),donors=';'.join(paired),eligibility='ELIGIBLE' if len(paired)>=3 else 'NE_FEWER_THAN_3_MATCHED_DONORS',action='REUSE_TARGETED' if (ligand,ct) in oldkeys else 'FIT_NEW' if len(paired)>=3 else 'NO_MODEL',dose='See preserved publisher condition metadata; not an efficacy-normalized dose comparison',source=str(E/'source_data/all_released_group_cell_counts.tsv')))
save(pd.DataFrame(eligibility),'inputs/all_ligand_pool_eligibility.tsv');comp=pd.DataFrame(comparisons);save(comp,'inputs/comparison_inventory.tsv')
print('INVENTORY',comp.groupby(['condition_class','action']).size().to_dict(),flush=True)

# Read already extracted/validated counts, not the raw serialized object and not
# a new cell-QC audit. Aggregate only new required donor/pool/family/arm groups.
with (RAW/'counts_parsed.pkl').open('rb') as f: parsed=pickle.load(f)
attrs=pairs(parsed.a)
def mmap_slot(n):
 v=attrs[n].v
 return np.memmap(v['file'],dtype=v['dtype'],mode='r',shape=(v['length'],)) if isinstance(v,dict) else v
ii=mmap_slot('i');xx=mmap_slot('x');ptr=np.load(RAW/'column_pointers.npy')
obs=pd.read_csv(RAW/'cell_metadata.tsv.gz',sep='\t',usecols=['cell_id','assignment','pool','celltype','perturbation'])
assert len(obs)+1==len(ptr)
obs['arm']=obs.perturbation.where(~obs.perturbation.isin(plan['control']),'Control')
aggregated={};used_positions=[]
for k,indices in obs.groupby(['assignment','pool','celltype','arm'],sort=True).indices.items():
 if k not in needed_groups:continue
 assert len(indices)==counts[k],k
 arr=np.zeros(len(genes),dtype=np.int64)
 for cell in indices:
  start,end=int(ptr[cell]),int(ptr[cell+1]);arr[np.asarray(ii[start:end],dtype=np.int64)]+=np.asarray(xx[start:end],dtype=np.int64)
 aggregated[k]=arr;used_positions.extend(indices.tolist())
 if len(aggregated)%100==0:print('NEW_GROUPS',len(aggregated),'/',len(needed_groups),flush=True)
assert set(aggregated)==needed_groups
samplemeta=pd.DataFrame(sample_specs);samplemeta.insert(0,'sample_id',[f'ALL_ORG_{i+1:04d}' for i in range(len(samplemeta))])
save(samplemeta,'inputs/all_ligand_pseudobulk_metadata.tsv')
all_scores=[];jobs=[]
for ligand,ct in itertools.product(ligands+['No_Cytomix'],families):
 b=samplemeta[(samplemeta.ligand==ligand)&(samplemeta.celltype==ct)].copy()
 if b.empty:continue
 if (ligand,ct) in oldkeys:
  sc=oldscores[(oldscores.ligand==ligand)&(oldscores.celltype==ct)].copy();sc['value_origin']='REUSED_TARGETED';all_scores.append(sc);continue
 cols=[]
 for r in b.itertuples():
  match=oldmeta[(oldmeta.ligand==ligand)&(oldmeta.celltype==ct)&(oldmeta.donor==r.donor)&(oldmeta.treatment==r.treatment)]
  if len(match):
   assert len(match)==1 and match.iloc[0].pools==r.pools
   cols.append(oldcounts[:,match.index[0]].toarray().ravel())
  else:cols.append(sum(aggregated[(r.donor,p,ct,'Control' if r.treatment=='Control' else ligand)] for p in r.pools.split(';')))
 cts=np.stack(cols,axis=1).astype(np.int64);assert cts.shape==(len(genes),len(b)) and (cts>=0).all()
 ranks=rankdata(cts,axis=0,method='average')/len(genes)
 for m,ix in sets.items():
  sc=b.copy();sc['module']=m;sc['raw_mean_rank']=ranks[ix,:].mean(axis=0);sc['value_origin']='NEW_SAME_SAVED_COUNT_RULE';sc['n_programme_members']=len(ix);all_scores.append(sc)
 if b.donor.nunique()>=3:
  stem=ligand+'__'+ct;df=pd.DataFrame(cts,columns=b.sample_id);df.insert(0,'gene',genes)
  dest=ROOT/'inputs/model_counts'/(stem+'.tsv.gz');df.to_csv(dest,sep='\t',index=False,compression={'method':'gzip','compresslevel':1})
  save(b,'inputs/model_counts/'+stem+'_metadata.tsv');jobs.append(dict(ligand=ligand,celltype=ct,stem=stem,counts=str(dest),metadata=str(ROOT/'inputs/model_counts'/(stem+'_metadata.tsv'))))
save(pd.concat(all_scores,ignore_index=True),'results/all_ligand_sample_scores.tsv');save(pd.DataFrame(jobs),'inputs/new_model_jobs.tsv')
register([(FROZEN,'Frozen programme members'),(E/'review/analysis_plan_locked.json','Approved original model and scoring specification'),(E/'source_data/all_released_group_cell_counts.tsv','Saved complete group coverage'),(E/'source_data/pseudobulk_counts.mtx','Existing paired pseudobulks reused'),(E/'source_data/pseudobulk_metadata.tsv','Existing matching keys'),(RAW/'counts_parsed.pkl','Saved sparse-count descriptors'),(RAW/'cell_metadata.tsv.gz','Saved released cell annotations'),(RAW/'column_pointers.npy','Saved sparse offsets'),(RAW/'genes.tsv','Saved gene order')]+[(Path(attrs[n].v['file']),'Saved immutable sparse counts') for n in ['i','x']]+[(SAVED/n,'Reused targeted result') for n in ['module_GSEA.tsv','donor_module_scores.tsv','donor_paired_effects.tsv','module_coverage.tsv','model_inventory.tsv']]+[(CODEBASE/'analysis/organoid/analyze_organoid.R','Approved model code reused'),(TASK/'csbj_20261016/Manuscript_CSBJ_revised_clean.docx','Current manuscript context')]+[(TASK/n,'User instructions') for n in ['Codex_Incremental_Analysis_and_Figures.md','Manuscript_Text_and_Figure4_Revision.md']])
jw(ROOT/'review/preparation_summary.json',dict(status='PASS',ligands=len(ligands),planned_comparisons=len(comp),eligible=int(comp.eligibility.eq('ELIGIBLE').sum()),new_models=len(jobs),reused_models=len(oldkeys),partial_observed_donors_retained=True,only_new_groups_aggregated=len(aggregated),new_unique_cells_read=len(set(used_positions)),completed_utc=now()))
print('PREPARATION_COMPLETE',len(jobs),'new models;',len(oldkeys),'reused models',flush=True)
