#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(cowplot);library(png)})
a<-commandArgs(TRUE);stopifnot(length(a)==1);o<-normalizePath(a[1],winslash='/');x<-file.path(o,'spatial_extension');sd<-file.path(x,'source_data')
rd<-function(n)read.delim(file.path(sd,n),check.names=FALSE,stringsAsFactors=FALSE)
pt<-rd('patient_scores.tsv');tt<-rd('frozen_transfer_results.tsv');ss<-rd('sample_manifest.tsv');qc<-rd('section_QC.tsv');sp<-rd('spot_scores.tsv.gz')
mods<-c('CURATED_IFN_VISIBILITY','ICI_COLITIS_EPITHELIAL_UP_STRICT','ICI_COLITIS_EPITHELIAL_UP_NON_IFN','ICI_COLITIS_EPITHELIAL_LOSS_STRICT')
labels<-c('IFN visibility','Disease-up','Non-IFN residual','Disease-depleted')
pt$condition<-factor(pt$condition,levels=c('Healthy','ICI colitis','UC'));pt$module<-factor(pt$module,levels=mods)
group_n<-sapply(levels(pt$condition),function(g)length(unique(pt$patient[pt$condition==g])))
group_labels<-setNames(paste0(names(group_n),'\nn=',group_n),names(group_n))
theme_set(theme_classic(base_size=9)+theme(text=element_text(colour='#22313D'),axis.text=element_text(colour='#22313D'),plot.title=element_text(face='bold',size=10),plot.subtitle=element_text(size=8),plot.margin=margin(7,7,7,7)))
dimensions<-list();maps<-list()
savep<-function(p,id,h=7.2){
 base<-file.path(x,'figures',id)
 ggsave(paste0(base,'.pdf'),p,width=7.2,height=h,device=cairo_pdf,bg='white')
 ggsave(paste0(base,'.tiff'),p,width=7.2,height=h,dpi=600,device=grDevices::tiff,compression='lzw',bg='white')
 ggsave(paste0(base,'_preview.png'),p,width=7.2,height=h,dpi=200,bg='white')
 dimensions[[length(dimensions)+1]]<<-data.frame(figure=id,width_in=7.2,height_in=h,tiff_dpi=600)
 message('Saved ',id)
}
plots<-lapply(seq_along(mods),function(i){
 d<-pt[pt$module==mods[i],];t<-tt[tt$module==mods[i],]
 ggplot(d,aes(condition,oriented_score,colour=condition))+stat_summary(fun=median,geom='crossbar',width=.5,linewidth=.35,show.legend=FALSE)+geom_point(position=position_jitter(width=.08,height=0,seed=20260826),size=2.4,show.legend=FALSE)+scale_colour_manual(values=c('#2C79A6','#C14A43','#959EA5'))+labs(title=labels[i],subtitle=sprintf('ICI vs healthy: P=%.3g; FDR=%.3g',t$P,t$FDR_new_four_module_family),x=NULL,y=if(i%%2)'Disease-oriented mean gene rank'else NULL)+scale_x_discrete(labels=group_labels)
})
panel<-plot_grid(plotlist=plots,ncol=2,labels=letters[1:4],label_size=12)
head<-ggdraw()+draw_label('Frozen programmes in an additional spatial cohort',x=.02,hjust=0,fontface='bold',size=12)
foot<-ggdraw()+draw_label(sprintf('One point per donor; repeated sections aggregated before scoring.\n%s sections from %s donors; UC provides context without a new hypothesis test.\nHigher oriented score = expected disease direction (disease-depleted: 1 - raw rank).\nWhole-tissue measurements; sex, source batch and composition remain potential confounders.',nrow(ss),length(unique(ss$patient))),x=.03,hjust=0,size=8)
savep(plot_grid(head,panel,foot,ncol=1,rel_heights=c(.32,5.8,.85)),'Supplementary_Figure_8_new_spatial_transfer',7.2)
write.table(pt,file.path(sd,'SF8_patient_points.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(tt,file.path(sd,'SF8_comparison_labels.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
showmods<-mods[c(1,3,4)];showlabs<-c('IFN visibility','Non-IFN residual','Disease-depleted')
ranges<-lapply(showmods,function(m)range(sp$raw_score[sp$module==m],finite=TRUE));names(ranges)<-showmods
palette<-c('#F4F4EE','#B4D7D2','#4FA39F','#236F83','#173F5F')
one<-function(s,m=NULL,title=''){
 q<-qc[qc$section==s,];row<-ss[ss$section==s,];im<-png::readPNG(row$image);w<-dim(im)[2];h<-dim(im)[1]
 p<-ggplot()+annotation_raster(im,xmin=0,xmax=w,ymin=-h,ymax=0,interpolate=TRUE)
 if(!is.null(m)){
  d<-sp[sp$section==s&sp$module==m,];stopifnot(nrow(d)==q$kept_tissue_spots,all(is.finite(d$x)),all(is.finite(d$y)))
  p<-p+geom_point(data=d,aes(x,-y,colour=raw_score),size=.65,alpha=.9)+scale_colour_gradientn(colours=palette,limits=ranges[[m]],name='Raw mean rank',na.value='#999999')
 }
 p+coord_fixed(xlim=c(0,w),ylim=c(-h,0),expand=FALSE)+theme_void()+theme(plot.title=element_text(size=8,face='bold',hjust=.5),legend.position='none',plot.margin=margin(2,2,2,2))+labs(title=title)
}
legendplots<-lapply(seq_along(showmods),function(i){
 d<-data.frame(x=c(0,1),z=ranges[[showmods[i]]]);p<-ggplot(d,aes(x,x,colour=z))+geom_point()+scale_colour_gradientn(colours=palette,limits=ranges[[showmods[i]]],name=paste0(showlabs[i],' rank'),breaks=ranges[[showmods[i]]],labels=sprintf('%.2f',ranges[[showmods[i]]]))+theme_void()+theme(legend.position='bottom',legend.title=element_text(size=7),legend.text=element_text(size=7))+guides(colour=guide_colourbar(title.position='top',barwidth=unit(2.5,'cm'),barheight=unit(.2,'cm')))
 get_legend(p)
})
for(k in seq_len(ceiling(nrow(ss)/3))){
 chunk<-ss[seq((k-1)*3+1,min(k*3,nrow(ss))),];panels<-list()
 for(j in seq_len(nrow(chunk))){
  r<-chunk[j,];panels[[length(panels)+1]]<-one(r$section,title=paste0(r$section,' | ',r$condition,'\n',r$patient,' | H&E'))
  for(i in seq_along(showmods))panels[[length(panels)+1]]<-one(r$section,showmods[i],if(j==1)paste0(showlabs[i],'\nexpression rank')else'')
 }
 title<-ggdraw()+draw_label(paste0('Spatial distribution of frozen programmes | plate ',LETTERS[k]),x=.015,hjust=0,fontface='bold',size=12)
 body<-plot_grid(plotlist=panels,ncol=4,align='hv')
 legends<-plot_grid(NULL,legendplots[[1]],legendplots[[2]],legendplots[[3]],ncol=4)
 caption<-ggdraw()+draw_label('Registered native tissue images and released coordinates; all eligible sections are shown.\nScales are fixed across all sections within each module; rows do not share a physical zoom.\nHigher raw rank = higher expression. Spots are not independent patients; no spatial P values.',x=.02,hjust=0,size=7.5)
 id<-paste0('Supplementary_Figure_9',LETTERS[k],'_spatial_atlas')
 savep(plot_grid(title,body,legends,caption,ncol=1,rel_heights=c(.4,5.45,.45,.6)),id,7.1)
 maps[[k]]<-data.frame(figure=id,sections=paste(chunk$section,collapse=';'),sources='spot_scores.tsv.gz; section_QC.tsv; sample_manifest.tsv',record_key='GSM accession + barcode + module',unit='raw within-spot gene expression rank; image pixel coordinates')
}
write.table(do.call(rbind,dimensions),file.path(x,'review/figure_dimensions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
maps[[length(maps)+1]]<-data.frame(figure='Supplementary_Figure_8_new_spatial_transfer',sections='panels a-d; all eligible patients',sources='patient_scores.tsv; frozen_transfer_results.tsv; module_coverage.tsv',record_key='patient + module; test record_key',unit='patient after pooling repeated sections')
write.table(do.call(rbind,maps),file.path(x,'review/panel_source_map.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(x,'review/R_sessionInfo.txt'))
