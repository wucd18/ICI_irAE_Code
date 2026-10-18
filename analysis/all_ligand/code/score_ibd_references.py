from common import *
from scipy.io import mmread
from scipy.stats import rankdata
import shutil
cfg=jr(ROOT/'analysis_config.json')
defs=pd.read_csv(ROOT/'inputs/IBD_membership_definitions.tsv',sep='\t')
genes=pd.read_csv(ROOT/'inputs/pseudobulk_genes.tsv',sep='\t').gene.tolist();gm={g:i for i,g in enumerate(genes)}
sets={m:[gm[g] for g in b.gene if g in gm] for m,b in defs.groupby('module')}
coverage=[dict(module=m,requested=len(b),raw_available=len(sets[m]),missing_members=';'.join(g for g in b.gene if g not in gm),scoring='Mean within-sample expression ranks across 36912 genes; no control-gene subtraction, no learned weights',reference_type='Published membership set, not reproduction of publisher module scores') for m,b in defs.groupby('module')]
save(pd.DataFrame(coverage),'results/IBD_reference_coverage.tsv')
jw(ROOT/'inputs/IBD_reference_config.json',dict(locked_utc=now(),definition_source='Published Supplementary Data 2, sheet B Genes in IBD1 and IBD2',source_archive=str(Path(r'E:\ICI_irAE_Original_Project\01_data\GSE313368_article_supplementary.zip')),gene_set_semantics='Published member lists are used as sets; CA1 repeated twice in IBD2 is counted once. These are equal-weight membership benchmarks, explicitly not reconstructed cNMF or publisher Module_C1/Module_C2 activity scores.',direction='Raw higher expression rank; no sign inversion of IBD1 or IBD2',weights='All unique members weight 1, as membership sets only',symbol_matching='Exact symbols; missing symbols reported; no alias remapping',benchmark='Eligible matched donor pseudobulks only; same raw-rank scoring as ICI/IFN readouts',source_code_evidence='https://raw.githubusercontent.com/Genentech/secretome_dictionary/3c80f66552d9c3ea344ba02304eb588b717e9ed8/Capeling_Figure_1_Code.Rmd lines 147-149 identify IBD1/IBD2 member lists; original module scores at lines 289-306 are a different method'))
scores=[]
def append_scores(cts,b):
    assert cts.shape[0]==len(genes)
    ranks=rankdata(cts,axis=0,method='average')/len(genes)
    for m,ix in sets.items():
        s=b.copy();s['module']=m;s['raw_mean_rank']=ranks[ix,:].mean(axis=0);s['n_programme_members']=len(ix);s['value_origin']='NEW_IBD_MEMBERSHIP_RANK_ON_SAME_INPUT';scores.append(s)
oldm=pd.read_csv(E/'source_data/pseudobulk_metadata.tsv',sep='\t');oldc=mmread(E/'source_data/pseudobulk_counts.mtx').tocsc()
inv=pd.read_csv(ROOT/'inputs/comparison_inventory.tsv',sep='\t')
for r in inv[inv.action.eq('REUSE_TARGETED')].itertuples():
    b=oldm[(oldm.ligand==r.ligand)&(oldm.celltype==r.celltype)].copy();append_scores(oldc[:,b.index].toarray(),b)
jobs=pd.read_csv(ROOT/'inputs/new_model_jobs.tsv',sep='\t')
for i,r in enumerate(jobs.itertuples()):
    b=pd.read_csv(r.metadata,sep='\t');c=pd.read_csv(r.counts,sep='\t');assert c.gene.tolist()==genes and c.columns[1:].tolist()==b.sample_id.tolist();append_scores(c.iloc[:,1:].to_numpy(),b)
    if i%50==0:print('IBD_MODELS_SCORED',i,flush=True)
save(pd.concat(scores,ignore_index=True),'results/IBD_reference_sample_scores.tsv')
register([(Path(r'E:\ICI_irAE_Original_Project\01_data\GSE313368_article_supplementary.zip'),'Publisher IBD member sets and verbatim protein doses')])
print('IBD_COMPLETE',coverage,flush=True)
