from common import *
import shutil
dest=ROOT.with_name(ROOT.name.replace('_r2','_r3'));assert not dest.exists()
record=dict(status='FAILED_ENRICHMENT_PARSER',cause='UTF-8 BOM caused R read.csv under LC_ALL=C to call the first header X...module; all requested gene sets were empty. This is a technical failure, not biological NE.',evidence='read.csv names X...module; read.csv(fileEncoding=UTF-8-BOM) names module',fitted_gene_models=257,fit_results_reusable='All 257 gene-statistic files and model inventories are intact; same approved DESeq2 code; validate before use',next_run=str(dest),action='Reuse the new fits and donor scores, repair only enrichment from preserved ranks; no model refitting',R_exit_issue='Every R worker emitted SHARD_COMPLETE then returned Windows 0xC0000005; recorded separately, not marked process PASS',created=now())
jw(ROOT/'review/RECOVERY_REQUIRED.json',record)
dest.mkdir()
for sub in ['code','inputs','text','review','figure_source_data']:
    shutil.copytree(ROOT/sub,dest/sub,ignore=shutil.ignore_patterns('__pycache__'))
(dest/'logs').mkdir();shutil.copytree(ROOT/'logs',dest/'logs/reused_model_fit_logs')
(dest/'results/new_models').mkdir(parents=True);(dest/'figures').mkdir();(dest/'temp').mkdir()
for p in (ROOT/'results').iterdir():
    if p.is_file():shutil.copy2(p,dest/'results'/p.name)
for p in (ROOT/'results/new_models').iterdir():
    if p.name.endswith(('_genes.tsv.gz','_model.tsv')):shutil.copy2(p,dest/'results/new_models'/p.name)
    elif p.name.endswith('_complete.json'):shutil.copy2(p,dest/'results/new_models'/p.name.replace('_complete.json','_fit_reused.json'))
shutil.copy2(ROOT/'analysis_config.json',dest/'analysis_config.json')
jobs=pd.read_csv(dest/'inputs/new_model_jobs.tsv',sep='\t')
for c in ['counts','metadata']:jobs[c]=jobs[c].str.replace(str(ROOT),str(dest),regex=False)
jobs.to_csv(dest/'inputs/new_model_jobs.tsv',sep='\t',index=False)
for i in range(4):jobs.iloc[i::4].to_csv(dest/'inputs'/f'jobs_{i}.tsv',sep='\t',index=False)
code=dest/'code/increment_models.R';text=code.read_text();text=text.replace("read.csv(file.path(O,'inputs/frozen_module_definitions.csv'))","read.csv(file.path(O,'inputs/frozen_module_definitions.csv'),fileEncoding='UTF-8-BOM')\nstopifnot(all(c('module','gene')%in%names(defs)))")
text=text.replace("names(sets)<-cfg$modules_in_original_fgsea_order","names(sets)<-cfg$modules_in_original_fgsea_order\nstopifnot(all(lengths(sets)>0))")
code.write_text(text,encoding='utf8')
jw(dest/'review/recovery_provenance.json',record)
print('RECOVERY_CREATED',dest,flush=True)
