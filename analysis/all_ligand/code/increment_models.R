suppressPackageStartupMessages({library(DESeq2);library(fgsea);library(jsonlite)})
options(stringsAsFactors=FALSE,warn=1)
args<-commandArgs(TRUE);stopifnot(length(args)==2)
O<-normalizePath(args[1],winslash='/');shard<-args[2]
cfg<-fromJSON(file.path(O,'analysis_config.json'))
fit_paired<-function(counts,meta) {
 stopifnot(!anyDuplicated(paste(meta$donor,meta$treatment)),all(table(meta$donor)==2))
 meta$donor<-factor(meta$donor);meta$treatment<-factor(meta$treatment,levels=c('Control','Treatment'))
 design<-model.matrix(~donor+treatment,meta);stopifnot(qr(design)$rank==ncol(design))
 keep<-rowSums(counts)>=10
 dds<-DESeqDataSetFromMatrix(counts[keep,,drop=FALSE],meta,design=~donor+treatment)
 dds<-DESeq(dds,quiet=TRUE,parallel=FALSE,betaPrior=FALSE,minReplicatesForReplace=Inf)
 res<-results(dds,contrast=c('treatment','Treatment','Control'),alpha=.05,independentFiltering=FALSE)
 list(dds=dds,result=as.data.frame(res),keep=keep,design=design)
}
if(shard=='fixture') {
 set.seed(20260826)
 meta<-data.frame(donor=rep(c('a','b','c'),each=2),treatment=rep(c('Control','Treatment'),3),row.names=paste0('s',1:6))
 means<-matrix(rep(exp(seq(log(25),log(200),length.out=600)),6),600,6);means[1:75,c(2,4,6)]<-means[1:75,c(2,4,6)]*3
 counts<-matrix(rnbinom(length(means),mu=means,size=12),nrow=600,dimnames=list(paste0('g',1:600),rownames(meta)))
 a<-fit_paired(counts,meta)
 stopifnot(median(a$result$log2FoldChange[1:75])>1,abs(median(a$result$log2FoldChange[76:600]))<.3)
 ranks<-a$result$stat;names(ranks)<-rownames(a$result);set.seed(20260826)
 g<-fgseaSimple(list(fixture=paste0('g',1:75)),sort(ranks,decreasing=TRUE),nperm=1000,minSize=10,maxSize=500,nproc=1)
 stopifnot(g$NES>0,is.finite(g$pval))
 write_json(list(status='PASS',scope='synthetic paired model and enrichment, no old biological result reproduced'),file.path(O,'review/R_fixture.json'),auto_unbox=TRUE,pretty=TRUE)
 cat('FIXTURE_COMPLETE\n');quit(status=0)
}
jobs<-read.delim(file.path(O,'inputs',paste0('jobs_',shard,'.tsv')),check.names=FALSE)
defs<-read.csv(file.path(O,'inputs/frozen_module_definitions.csv'),fileEncoding='UTF-8-BOM')
stopifnot(all(c('module','gene')%in%names(defs)))
sets<-lapply(cfg$modules_in_original_fgsea_order,function(m)unique(defs$gene[defs$module==m]));names(sets)<-cfg$modules_in_original_fgsea_order
stopifnot(all(lengths(sets)>0))
capture.output(sessionInfo(),file=file.path(O,'logs',paste0('R_session_',shard,'.txt')))
for(j in seq_len(nrow(jobs))) {
 job<-jobs[j,];dest<-file.path(O,'results/new_models',job$stem)
 if(file.exists(paste0(dest,'_complete.json')))stop('Refusing to overwrite finished model')
 started<-format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z')
 cat('MODEL_START',job$stem,'\n');flush.console()
 b<-read.delim(job$metadata,check.names=FALSE);d<-read.delim(gzfile(job$counts),check.names=FALSE)
 genes<-d$gene;cts<-as.matrix(d[,-1,drop=FALSE]);storage.mode(cts)<-'integer';rownames(cts)<-genes;rownames(b)<-b$sample_id
 stopifnot(identical(colnames(cts),b$sample_id),length(unique(b$donor))==3)
 a<-fit_paired(cts,b);res<-a$result;res$gene<-rownames(res)
 full<-data.frame(gene=genes,raw_total=rowSums(cts),retained_total_count_ge10=a$keep)
 full<-merge(full,res,by='gene',all.x=TRUE,sort=FALSE)
 full$maxCooks<-mcols(a$dds)$maxCooks[match(full$gene,rownames(a$dds))]
 full$rank_eligible<-is.finite(full$stat)&is.finite(full$pvalue)
 full$ligand<-job$ligand;full$celltype<-job$celltype;full$n_donors<-length(unique(b$donor))
 z<-full[full$rank_eligible,];ranks<-z$stat;names(ranks)<-z$gene;ranks<-sort(ranks,decreasing=TRUE)
 stopifnot(length(ranks)>=1000,!anyDuplicated(names(ranks)))
 set.seed(20260826)
 fg<-as.data.frame(fgseaSimple(sets,ranks,nperm=10000,minSize=10,maxSize=500,scoreType='std',nproc=1))
 records<-list();coverage<-list()
 for(m in names(sets)) {
  v<-intersect(sets[[m]],names(ranks));f<-fg[fg$pathway==m,]
  coverage[[length(coverage)+1]]<-data.frame(ligand=job$ligand,celltype=job$celltype,module=m,requested=length(sets[[m]]),raw_available=length(intersect(sets[[m]],genes)),rank_available=length(v))
  records[[length(records)+1]]<-data.frame(ligand=job$ligand,celltype=job$celltype,module=m,n_donors=length(unique(b$donor)),status=if(nrow(f)==1)'ESTIMATED' else 'NOT_ESTIMABLE_MODULE_COVERAGE',NES=if(nrow(f)==1)f$NES else NA_real_,pval=if(nrow(f)==1)f$pval else NA_real_,ES=if(nrow(f)==1)f$ES else NA_real_,size=length(v),leadingEdge=if(nrow(f)==1)paste(f$leadingEdge[[1]],collapse=';') else NA_character_,n_ranked_genes=length(ranks))
 }
 out<-do.call(rbind,records)
 stopifnot(nrow(out)==length(sets),all(is.finite(out$NES[out$status=='ESTIMATED'])))
 con<-gzfile(paste0(dest,'_genes.tsv.gz'),'wt');write.table(full,con,sep='\t',row.names=FALSE,quote=FALSE,na='NA');close(con)
 write.table(out,paste0(dest,'_GSEA.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
 write.table(do.call(rbind,coverage),paste0(dest,'_coverage.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
 sf<-sizeFactors(a$dds)
 inventory<-data.frame(ligand=job$ligand,celltype=job$celltype,n_donors=length(unique(b$donor)),n_pseudobulks=nrow(b),design_columns=paste(colnames(a$design),collapse=';'),n_filtered_genes=sum(a$keep),n_ranked_genes=length(ranks),n_cooks_or_P_missing=sum(is.finite(full$stat)&!is.finite(full$pvalue)),dispersion_fit_type=attr(dispersionFunction(a$dds),'fitType'),min_size_factor=min(sf),max_size_factor=max(sf))
 write.table(inventory,paste0(dest,'_model.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
 write_json(list(status='COMPLETE',model=job$stem,started=started,finished=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z')),paste0(dest,'_complete.json'),pretty=TRUE,auto_unbox=TRUE)
 cat('MODEL_COMPLETE',job$stem,'\n');flush.console()
 rm(a,res,full,cts,d);gc(verbose=FALSE)
}
write_json(list(status='COMPLETE',shard=shard,models=nrow(jobs)),file.path(O,'logs',paste0('shard_',shard,'_complete.json')),pretty=TRUE,auto_unbox=TRUE)
cat('SHARD_COMPLETE',shard,'\n')
