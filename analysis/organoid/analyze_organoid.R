suppressPackageStartupMessages({library(Matrix);library(DESeq2);library(fgsea);library(jsonlite)})
args<-commandArgs(TRUE);stopifnot(length(args)==1);O<-normalizePath(args[1],winslash='/');E<-file.path(O,'organoid_donor_extension')
options(stringsAsFactors=FALSE,warn=1)
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
if(identical(Sys.getenv('ICI_FIXTURE_ONLY'),'1')) {
 set.seed(20260826);meta<-data.frame(donor=rep(c('a','b','c'),each=2),treatment=rep(c('Control','Treatment'),3),row.names=paste0('s',1:6))
 means<-matrix(rep(exp(seq(log(25),log(200),length.out=600)),6),600,6)
 means[1:75,c(2,4,6)]<-means[1:75,c(2,4,6)]*3
 counts<-matrix(rnbinom(length(means),mu=means,size=12),nrow=600,dimnames=list(paste0('g',1:600),rownames(meta)))
 a<-fit_paired(counts,meta)
 stopifnot(median(a$result$log2FoldChange[1:75])>1,abs(median(a$result$log2FoldChange[76:600]))<.3,identical(colnames(a$design),c('(Intercept)','donorb','donorc','treatmentTreatment')))
 ranks<-a$result$stat;names(ranks)<-rownames(a$result);set.seed(20260826);g<-fgseaSimple(list(fixture=paste0('g',1:75)),sort(ranks,decreasing=TRUE),nperm=1000,minSize=10,maxSize=500,nproc=1)
 stopifnot(g$NES>0,is.finite(g$pval));cat('PAIRED_DESEQ_GSEA_FIXTURE_PASS\n');quit(status=0)
}
plan<-fromJSON(file.path(E,'review/analysis_plan_locked.json'));stopifnot(fromJSON(file.path(E,'review/preflight_QC.json'))$status=='PASS')
meta<-read.delim(file.path(E,'source_data/pseudobulk_metadata.tsv'),check.names=FALSE)
genes<-read.delim(file.path(E,'source_data/pseudobulk_genes.tsv'))$gene
counts<-readMM(file.path(E,'source_data/pseudobulk_counts.mtx'));rownames(counts)<-genes;colnames(counts)<-meta$sample_id
defs<-read.delim(file.path(O,'supplementary/original_tables/Table_S21_frozen_module_definitions.tsv'))
sets<-lapply(plan$modules,function(m)unique(defs$gene[defs$module==m]));names(sets)<-plan$modules
outdir<-file.path(E,'results',format(Sys.time(),'%Y%m%dT%H%M%S'));dir.create(outdir,recursive=TRUE);gene_dir<-file.path(outdir,'gene_results');dir.create(gene_dir)
summaries<-list();score_rows<-list();score_tests<-list();coverage<-list();models<-list()
reference<-trimws(strsplit(plan$assay_reference,';',fixed=TRUE)[[1]][1])
conditions<-c(plan$primary_ligand,plan$secondary_ligands,reference)
stopifnot(all(conditions %in% meta$ligand),!anyDuplicated(conditions))
for(ligand in conditions)for(ct in plan$celltypes) {
 cat('MODEL',ligand,ct,'\n');flush.console()
 b<-meta[meta$ligand==ligand&meta$celltype==ct,];n<-length(unique(b$donor));family<-if(ligand==plan$primary_ligand)'primary_TNFSF12' else if(ligand==reference)'assay_reference' else 'secondary_original_nominees'
 if(n<3) {
  for(m in names(sets))summaries[[length(summaries)+1]]<-data.frame(ligand=ligand,celltype=ct,module=m,family=family,n_donors=n,status='NOT_ESTIMABLE_FEWER_THAN_3_DONORS',NES=NA_real_,pval=NA_real_,ES=NA_real_,size=NA_integer_,leadingEdge=NA_character_,n_ranked_genes=NA_integer_)
  next
 }
 cts<-as.matrix(counts[,b$sample_id,drop=FALSE]);storage.mode(cts)<-'integer';rownames(b)<-b$sample_id
 a<-fit_paired(cts,b);res<-a$result;res$gene<-rownames(res)
 full<-data.frame(gene=genes,raw_total=rowSums(cts),retained_total_count_ge10=a$keep)
 full<-merge(full,res,by='gene',all.x=TRUE,sort=FALSE)
 full$maxCooks<-mcols(a$dds)$maxCooks[match(full$gene,rownames(a$dds))]
 full$rank_eligible<-is.finite(full$stat)&is.finite(full$pvalue)
 full$ligand<-ligand;full$celltype<-ct;full$n_donors<-n
 con<-gzfile(file.path(gene_dir,paste0(ligand,'__',ct,'.tsv.gz')),'wt');write.table(full,con,sep='\t',row.names=FALSE,quote=FALSE,na='NA');close(con)
 z<-full[full$rank_eligible,];ranks<-z$stat;names(ranks)<-z$gene;ranks<-sort(ranks,decreasing=TRUE)
 stopifnot(length(ranks)>=1000,!anyDuplicated(names(ranks)))
 set.seed(20260826)
 fg<-as.data.frame(fgseaSimple(sets,ranks,nperm=10000,minSize=10,maxSize=500,scoreType='std',nproc=1))
 for(m in names(sets)) {
  cov<-intersect(sets[[m]],genes);v<-intersect(sets[[m]],names(ranks));coverage[[length(coverage)+1]]<-data.frame(ligand=ligand,celltype=ct,module=m,requested=length(sets[[m]]),raw_available=length(cov),rank_available=length(v))
  f<-fg[fg$pathway==m,]
  if(nrow(f)==1)summaries[[length(summaries)+1]]<-data.frame(ligand=ligand,celltype=ct,module=m,family=family,n_donors=n,status='ESTIMATED',NES=f$NES,pval=f$pval,ES=f$ES,size=f$size,leadingEdge=paste(f$leadingEdge[[1]],collapse=';'),n_ranked_genes=length(ranks))
  else summaries[[length(summaries)+1]]<-data.frame(ligand=ligand,celltype=ct,module=m,family=family,n_donors=n,status='NOT_ESTIMABLE_MODULE_COVERAGE',NES=NA_real_,pval=NA_real_,ES=NA_real_,size=length(v),leadingEdge=NA_character_,n_ranked_genes=length(ranks))
  raw<-colMeans(apply(cts,2,rank,ties.method='average')[match(cov,genes),,drop=FALSE]/length(genes))
  sd<-data.frame(ligand=ligand,celltype=ct,module=m,sample_id=b$sample_id,donor=b$donor,treatment=b$treatment,raw_mean_rank=raw,n_cells=b$n_cells,pools=b$pools);score_rows[[length(score_rows)+1]]<-sd
  donors<-sort(unique(sd$donor));change<-sapply(donors,function(d)sd$raw_mean_rank[sd$donor==d&sd$treatment=='Treatment']-sd$raw_mean_rank[sd$donor==d&sd$treatment=='Control'])
  signs<-as.matrix(expand.grid(rep(list(c(-1,1)),length(change))));null<-as.vector(signs%*%change)/length(change);exact<-mean(abs(null)>=abs(mean(change))-1e-15)
  score_tests[[length(score_tests)+1]]<-data.frame(ligand=ligand,celltype=ct,module=m,family=family,n_donors=n,mean_paired_rank_change=mean(change),median_paired_rank_change=median(change),min_paired_change=min(change),max_paired_change=max(change),n_positive=sum(change>0),n_negative=sum(change<0),exact_signflip_P=exact)
 }
 sf<-sizeFactors(a$dds)
 models[[length(models)+1]]<-data.frame(ligand=ligand,celltype=ct,n_donors=n,n_pseudobulks=nrow(b),design_columns=paste(colnames(a$design),collapse=';'),n_filtered_genes=sum(a$keep),n_ranked_genes=length(ranks),n_cooks_or_P_missing=sum(is.finite(full$stat)&!is.finite(full$pvalue)),dispersion_fit_type=attr(dispersionFunction(a$dds),'fitType'),min_size_factor=min(sf),max_size_factor=max(sf))
 write_json(list(last_completed=paste(ligand,ct),models_completed=length(models),output=outdir),file.path(outdir,'analysis_progress.json'),pretty=TRUE,auto_unbox=TRUE)
 rm(a,full,res,cts);gc(verbose=FALSE)
}
r<-do.call(rbind,summaries);sc<-do.call(rbind,score_rows);st<-do.call(rbind,score_tests)
family_n<-c(primary_TNFSF12=length(plan$celltypes)*length(plan$modules),secondary_original_nominees=length(plan$secondary_ligands)*length(plan$celltypes)*length(plan$modules),assay_reference=length(plan$celltypes)*length(plan$modules))
r$FDR_predeclared_family<-NA_real_;st$FDR_predeclared_family<-NA_real_
for(f in names(family_n)){
 ix<-r$family==f;r$FDR_predeclared_family[ix]<-p.adjust(r$pval[ix],method='BH',n=family_n[[f]])
 iy<-st$family==f;st$FDR_predeclared_family[iy]<-p.adjust(st$exact_signflip_P[iy],method='BH',n=family_n[[f]])
}
save_t<-function(x,n)write.table(x,file.path(outdir,n),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
save_t(r,'module_GSEA.tsv');save_t(sc,'donor_module_scores.tsv');save_t(st,'donor_paired_effects.tsv');save_t(do.call(rbind,coverage),'module_coverage.tsv');save_t(do.call(rbind,models),'model_inventory.tsv')
capture.output(sessionInfo(),file=file.path(outdir,'sessionInfo.txt'))
write_json(list(status='PASS',output=outdir,models_completed=length(models),planned_GSEA_tests=nrow(r),estimated_GSEA_tests=sum(r$status=='ESTIMATED'),primary_GSEA_FDR05=sum(r$family=='primary_TNFSF12'&r$FDR_predeclared_family<.05,na.rm=TRUE),primary_exact_donor_FDR05=sum(st$family=='primary_TNFSF12'&st$FDR_predeclared_family<.05,na.rm=TRUE)),file.path(outdir,'analysis_complete.json'),pretty=TRUE,auto_unbox=TRUE)
print(r[r$family=='primary_TNFSF12',c('celltype','module','NES','pval','FDR_predeclared_family')]);cat('ORGANOID_DONOR_ANALYSIS_PASS\n')
