suppressPackageStartupMessages({library(fgsea);library(jsonlite)})
options(stringsAsFactors=FALSE,warn=1)
args<-commandArgs(TRUE);O<-normalizePath(args[1],winslash='/');shard<-args[2]
cfg<-fromJSON(file.path(O,'analysis_config.json'))
defs<-read.csv(file.path(O,'inputs/frozen_module_definitions.csv'),fileEncoding='UTF-8-BOM')
stopifnot(all(c('module','gene')%in%names(defs)))
sets<-lapply(cfg$modules_in_original_fgsea_order,function(m)unique(defs$gene[defs$module==m]));names(sets)<-cfg$modules_in_original_fgsea_order
expected<-fromJSON(file.path(O,'inputs/expected_module_member_counts.json'))
stopifnot(all(lengths(sets)>0),identical(as.integer(lengths(sets)),as.integer(unlist(expected[names(sets)]))))
capture.output(sessionInfo(),file=file.path(O,'logs',paste0('R_enrichment_session_',shard,'.txt')))
write.table(defs[defs$module%in%names(sets),c('module','gene')],file.path(O,'review',paste0('R_loaded_programme_members_',shard,'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
jobs<-read.delim(file.path(O,'inputs',paste0('jobs_',shard,'.tsv')),check.names=FALSE)
for(j in seq_len(nrow(jobs))) {
 job<-jobs[j,];dest<-file.path(O,'results/new_models',job$stem)
 if(file.exists(paste0(dest,'_complete.json')))stop('Refusing to overwrite completed enrichment')
 started<-format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z');cat('ENRICHMENT_START',job$stem,'\n');flush.console()
 full<-read.delim(gzfile(paste0(dest,'_genes.tsv.gz')),check.names=FALSE)
 stopifnot(!anyDuplicated(full$gene),all(full$rank_eligible==(is.finite(full$stat)&is.finite(full$pvalue))))
 z<-full[full$rank_eligible,];ranks<-z$stat;names(ranks)<-z$gene;ranks<-sort(ranks,decreasing=TRUE)
 stopifnot(length(ranks)>=1000,!anyDuplicated(names(ranks)))
 set.seed(20260826)
 fg<-as.data.frame(fgseaSimple(sets,ranks,nperm=10000,minSize=10,maxSize=500,scoreType='std',nproc=1))
 records<-list();coverage<-list()
 for(m in names(sets)) {
  v<-intersect(sets[[m]],names(ranks));f<-fg[fg$pathway==m,]
  stopifnot(length(v)>=10,nrow(f)==1,is.finite(f$NES),is.finite(f$pval))
  coverage[[length(coverage)+1]]<-data.frame(ligand=job$ligand,celltype=job$celltype,module=m,requested=length(sets[[m]]),raw_available=length(intersect(sets[[m]],full$gene)),rank_available=length(v))
  records[[length(records)+1]]<-data.frame(ligand=job$ligand,celltype=job$celltype,module=m,n_donors=length(unique(full$n_donors))*full$n_donors[1],status='ESTIMATED',NES=f$NES,pval=f$pval,ES=f$ES,size=length(v),leadingEdge=paste(f$leadingEdge[[1]],collapse=';'),n_ranked_genes=length(ranks))
 }
 out<-do.call(rbind,records);stopifnot(nrow(out)==4,all(out$n_donors==3))
 write.table(out,paste0(dest,'_GSEA.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
 write.table(do.call(rbind,coverage),paste0(dest,'_coverage.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
 fit<-fromJSON(paste0(dest,'_fit_reused.json'))
 write_json(list(status='COMPLETE',model=job$stem,started=fit$started,finished=fit$finished,enrichment_started=started,enrichment_finished=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'),fit_reuse='Unchanged new-model gene statistics from r2; no refit',parser='UTF-8-BOM explicitly decoded and gene sets cross-checked'),paste0(dest,'_complete.json'),pretty=TRUE,auto_unbox=TRUE)
 cat('ENRICHMENT_COMPLETE',job$stem,'\n');flush.console()
}
write_json(list(status='COMPLETE',shard=shard,models=nrow(jobs)),file.path(O,'logs',paste0('enrichment_',shard,'_complete.json')),pretty=TRUE,auto_unbox=TRUE)
cat('SHARD_COMPLETE',shard,'\n')
