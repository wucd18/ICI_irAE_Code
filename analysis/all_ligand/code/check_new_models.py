from common import *
from analysis_functions import *
import time
jobs=pd.read_csv(ROOT/'inputs/new_model_jobs.tsv',sep='\t');rows=[]
for i,r in enumerate(jobs.itertuples()):
    p=ROOT/'results/new_models'/r.stem
    # Wait only for genuinely running new models; no refitting or old-model audit.
    while not Path(str(p)+'_complete.json').exists():time.sleep(2)
    b=pd.read_csv(str(p)+'_genes.tsv.gz',sep='\t');mi=pd.read_csv(str(p)+'_model.tsv',sep='\t').iloc[0];gs=pd.read_csv(str(p)+'_GSEA.tsv',sep='\t');c=pd.read_csv(str(p)+'_coverage.tsv',sep='\t');meta=pd.read_csv(r.metadata,sep='\t')
    assert len(b)==36912 and b.gene.is_unique
    assert (b.raw_total>=0).all() and np.array_equal(b.retained_total_count_ge10.to_numpy(),(b.raw_total>=10).to_numpy())
    assert ((b.pvalue.dropna()>=0)&(b.pvalue.dropna()<=1)).all()
    assert np.array_equal(b.rank_eligible.to_numpy(),(np.isfinite(b.stat)&np.isfinite(b.pvalue)).to_numpy())
    assert b.rank_eligible.sum()==mi.n_ranked_genes and mi.n_donors==3 and mi.n_pseudobulks==6
    assert len(meta)==6 and meta.groupby('donor').size().eq(2).all() and meta.n_cells.ge(20).all()
    assert (c.rank_available<=c.raw_available).all() and (c.raw_available<=c.requested).all()
    assert np.array_equal(np.sign(gs.NES),np.sign(gs.ES))
    completion=jr(str(p)+'_complete.json')
    fit_root=ROOT.with_name(ROOT.name.replace('_r3','_r2'));old_gene=fit_root/'results/new_models'/(r.stem+'_genes.tsv.gz')
    assert sha(str(p)+'_genes.tsv.gz')==sha(old_gene)
    assert c.requested.min()>0 and c.rank_available.min()>=10
    rows.append(dict(ligand=r.ligand,celltype=r.celltype,status='PASS',genes=len(b),n_ranked_genes=mi.n_ranked_genes,n_pseudobulks=len(meta),count_input=r.counts,count_input_sha256=sha(r.counts),metadata_input=r.metadata,metadata_sha256=sha(r.metadata),fit_code_sha256=sha(fit_root/'code/increment_models.R'),enrichment_code_sha256=sha(ROOT/'code/recover_enrichment.R'),gene_fit_source=str(old_gene),gene_output_sha256=sha(str(p)+'_genes.tsv.gz'),GSEA_output_sha256=sha(str(p)+'_GSEA.tsv'),started=completion['started'],finished=completion['finished'],enrichment_started=completion['enrichment_started'],enrichment_finished=completion['enrichment_finished']))
    if (i+1)%50==0:print('NEW_MODEL_CHECKED',i+1,flush=True)
save(pd.DataFrame(rows),'review/new_model_checks_and_provenance.tsv')
defs=pd.read_csv(ROOT/'inputs/frozen_module_definitions.csv');cfg=jr(ROOT/'analysis_config.json');expected=set(map(tuple,defs[defs.module.isin(cfg['modules_in_original_fgsea_order'])][['module','gene']].to_numpy()))
for shard in range(4):assert set(map(tuple,pd.read_csv(ROOT/'review'/f'R_loaded_programme_members_{shard}.tsv',sep='\t').to_numpy()))==expected
print('NEW_MODEL_CHECKS_COMPLETE',len(rows),flush=True)
