suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(cowplot)})
a<-commandArgs(TRUE);stopifnot(length(a)==1);o<-normalizePath(a[1],winslash='/');x<-file.path(o,'xenium_extension');sd<-file.path(x,'source_data')
rd<-function(n)read.delim(file.path(sd,n),check.names=FALSE,stringsAsFactors=FALSE)
g<-rd('patient_gene_expression.tsv');p<-rd('patient_compartment_scores.tsv');cv<-rd('module_coverage.tsv');cells<-rd('cell_coordinates_scores.tsv.gz')
# Python emits True/False; R may retain these as character strings.
# Decode explicitly before subsetting so strings cannot become row names.
included<-match(tolower(as.character(cells$analysis_included)),c('false','true'))-1L
stopifnot(!anyNA(included),any(included==1L))
cells<-cells[included==1L,]
stopifnot(!anyNA(cells$cell_key),!anyDuplicated(cells$cell_key),all(is.finite(cells$x_centroid)),all(is.finite(cells$y_centroid)))
mods<-cv$module;modlab<-c('IFN visibility','Disease-up','Non-IFN residual','Disease-depleted');names(modlab)<-mods
comps<-c('IEC','Immune','Stromal','Endothelial','Enteric_glia');cpal<-setNames(c('#377CA5','#C3534E','#B49B5A','#528E82','#9F79AD'),comps)
cond<-c('HC','ICI','UC');colors<-setNames(c('#2C79A6','#C14A43','#929BA3'),cond)
g$compartment<-factor(g$compartment,levels=comps);p$compartment<-factor(p$compartment,levels=comps);g$condition<-factor(g$condition,levels=cond);p$condition<-factor(p$condition,levels=cond)
theme_set(theme_classic(base_size=9)+theme(text=element_text(colour='#22313D'),axis.text=element_text(colour='#22313D'),plot.title=element_text(face='bold',size=11),plot.subtitle=element_text(size=8),plot.margin=margin(7,7,7,7)))
dimensions<-list();maps<-list()
savep<-function(plot,id,h){b<-file.path(x,'figures',id);ggsave(paste0(b,'.pdf'),plot,width=7.2,height=h,device=cairo_pdf,bg='white');ggsave(paste0(b,'.tiff'),plot,width=7.2,height=h,dpi=600,device=grDevices::tiff,compression='lzw',bg='white');ggsave(paste0(b,'_preview.png'),plot,width=7.2,height=h,dpi=200,bg='white');dimensions[[length(dimensions)+1]]<<-data.frame(figure=id,width_in=7.2,height_in=h,tiff_dpi=600);message('Saved ',id)}
panels<-lapply(c('SORL1','TNFSF12','TNFRSF12A'),function(gene){d<-g[g$gene==gene,];ggplot(d,aes(compartment,mean_molecules_per_cell,colour=condition))+geom_point(position=position_jitter(width=.15,height=0,seed=20260826),size=2)+scale_colour_manual(values=colors,name='Condition')+scale_x_discrete(labels=c('Epithelium','Immune','Stromal','Endothelium','Enteric glia'))+labs(title=gene,x=NULL,y='Mean molecules / cell')})
head<-ggdraw()+draw_label('Nominated genes across tissue compartments',x=.02,hjust=0,fontface='bold',size=12)
foot<-ggdraw()+draw_label('One point per patient and published cell compartment. All eligible patients are shown.\nTargeted Xenium panel; source segmentation and cell-subtype composition affect these estimates.\nDonor IDs overlap the original discovery cohort; this is spatial follow-through, not independent validation.',x=.02,hjust=0,size=8)
savep(plot_grid(head,plot_grid(plotlist=panels,ncol=1,labels=letters[1:3]),foot,ncol=1,rel_heights=c(.35,6,.65)),'Supplementary_Figure_10_spatial_gene_compartments',7.2)
maps[[1]]<-data.frame(figure='Supplementary_Figure_10_spatial_gene_compartments',panel='a-c',source='patient_gene_expression.tsv',key='patient|compartment|gene',statistical_unit='patient; mean observed molecules per cell, no hypothesis test')
ep<-p[p$compartment=='IEC',];npat<-sapply(cond,function(cc)length(unique(ep$patient[ep$condition==cc])))
panels<-lapply(seq_along(mods),function(i){d<-ep[ep$module==mods[i],];coverage<-cv[cv$module==mods[i],];ggplot(d,aes(condition,oriented_rank,colour=condition))+geom_point(position=position_jitter(width=.08,height=0,seed=20260826),size=2.6,show.legend=FALSE)+stat_summary(fun=median,geom='crossbar',width=.4,linewidth=.35,show.legend=FALSE)+scale_colour_manual(values=colors)+scale_x_discrete(labels=setNames(paste0(cond,'\nn=',npat),cond))+labs(title=modlab[[mods[i]]],subtitle=sprintf('Targeted coverage: %s / %s frozen genes',coverage$available,coverage$requested),x=NULL,y='Disease-oriented panel rank')})
head<-ggdraw()+draw_label('Frozen programmes within annotated epithelium',x=.02,hjust=0,fontface='bold',size=12)
foot<-ggdraw()+draw_label('Counts pooled within each patient and the published IEC compartment before scoring.\nModule members are unchanged; only assayed genes can contribute. No missing gene is imputed.\nOverlapping discovery donors and partial targeted coverage preclude an independent full-module validation.',x=.02,hjust=0,size=8)
savep(plot_grid(head,plot_grid(plotlist=panels,ncol=2,labels=letters[1:4]),foot,ncol=1,rel_heights=c(.35,4.8,.65)),'Supplementary_Figure_11_epithelial_spatial_scores',6.1)
maps[[2]]<-data.frame(figure='Supplementary_Figure_11_epithelial_spatial_scores',panel='a-d',source='patient_compartment_scores.tsv; module_coverage.tsv',key='patient|IEC|module',statistical_unit='patient; within-pseudobulk rank on targeted panel, no hypothesis test')
showmods<-mods[c(3,4)];vals<-lapply(showmods,function(m)range(cells[[m]],finite=TRUE));names(vals)<-showmods
pal<-c('#F1F2EC','#B9DAD2','#69B2AA','#237B8A','#173F5F');smax<-max(log1p(cells$SORL1_molecules))
stopifnot(is.finite(smax),smax>0,all(vapply(vals,function(v)all(is.finite(v))&&diff(v)>0,logical(1))))
one<-function(d,type,title){
 base<-ggplot(d,aes(x_centroid,-y_centroid))
 if(type=='compartment')base<-base+geom_point(aes(colour=compartment),size=.08,stroke=0)+scale_colour_manual(values=cpal,drop=FALSE)
 else if(type=='SORL1')base<-base+geom_point(aes(colour=log1p(SORL1_molecules)),size=.1,stroke=0)+scale_colour_gradientn(colours=pal,limits=c(0,smax))
 else base<-base+geom_point(aes(colour=.data[[type]]),size=.1,stroke=0)+scale_colour_gradientn(colours=pal,limits=vals[[type]])
 base+coord_fixed()+theme_void()+theme(legend.position='none',plot.title=element_text(size=8,hjust=.5,face='bold'),plot.margin=margin(2,2,2,2))+labs(title=title)
}
legendplot<-function(lim,title){d<-data.frame(x=1:2,y=lim);get_legend(ggplot(d,aes(x,x,colour=y))+geom_point()+scale_colour_gradientn(colours=pal,limits=lim,breaks=lim,labels=sprintf('%.2f',lim),name=title)+theme_void()+theme(legend.position='bottom',legend.text=element_text(size=7),legend.title=element_text(size=7))+guides(colour=guide_colourbar(title.position='top',barwidth=unit(2.5,'cm'),barheight=unit(.2,'cm'))))}
lg<-list(legendplot(c(0,smax),'SORL1 log1p molecules'),legendplot(vals[[showmods[1]]],'Non-IFN raw rank'),legendplot(vals[[showmods[2]]],'Depleted raw rank'))
comp_legend<-get_legend(ggplot(data.frame(x=seq_along(comps),compartment=factor(comps,levels=comps)),aes(x,x,colour=compartment))+geom_point(size=2)+scale_colour_manual(values=cpal,name=NULL,labels=c('Epithelium','Immune','Stromal','Endothelium','Enteric glia'))+theme_void()+theme(legend.position='bottom',legend.text=element_text(size=7))+guides(colour=guide_legend(nrow=1)))
viewmaps<-list();viewkeys<-list()
for(k in seq_along(cond)){
 cc<-cond[k];patients<-sort(unique(cells$patient[cells$condition==cc]));pan<-list()
 for(patient in patients){
  d<-cells[cells$patient==patient,]
  # Layout only: remove the largest empty vertical strip when it exceeds
  # 10% of this patient's coordinate range. No expression enters this rule.
  yy<-sort(unique(d$y_centroid));gap<-diff(yy);cut<-if(max(gap)>.1*diff(range(yy)))mean(yy[which.max(gap)+0:1]) else Inf
  d$view_tile<-ifelse(d$y_centroid<=cut,1L,2L)
  for(tile in sort(unique(d$view_tile))){
   z<-d[d$view_tile==tile,];key<-paste(patient,tile,sep='|')
   viewmaps[[length(viewmaps)+1]]<-data.frame(view_key=key,patient=patient,tile=tile,n_cells=nrow(z),xmin=min(z$x_centroid),xmax=max(z$x_centroid),ymin=min(z$y_centroid),ymax=max(z$y_centroid),empty_strip_split=cut)
   viewkeys[[length(viewkeys)+1]]<-data.frame(cell_key=z$cell_key,patient=patient,view_tile=tile)
   pan<-c(pan,list(one(z,'compartment',paste0(patient,' | tile ',tile,'\nPublished compartments')),one(z,'SORL1','SORL1'),one(z,showmods[1],'Non-IFN raw rank'),one(z,showmods[2],'Disease-depleted raw rank')))
  }
 }
 header<-ggdraw()+draw_label(paste('Spatial follow-through |',cc,'| all assigned patients'),x=.015,hjust=0,fontface='bold',size=12)
 footer<-ggdraw()+draw_label('Coordinate tiles remove large empty strips; all assigned cells remain. Within-tile geometry is preserved.\nTiles use separate view extents and are not biological replicates; full bounds and cell keys are supplied.\nColour scales are shared across every patient. No spatial or disease-comparison P values.',x=.02,hjust=0,size=7.5)
 id<-paste0('Supplementary_Figure_12',LETTERS[k],'_spatial_cell_atlas')
 savep(plot_grid(header,plot_grid(plotlist=pan,ncol=4),comp_legend,plot_grid(NULL,lg[[1]],lg[[2]],lg[[3]],ncol=4),footer,ncol=1,rel_heights=c(.4,7,.35,.45,.6)),id,8.8)
 maps[[length(maps)+1]]<-data.frame(figure=id,panel=paste(patients,collapse=';'),source='cell_coordinates_scores.tsv.gz',key='cell_key; patient; source compartment',statistical_unit='released cell centroids and observed molecules or raw ranks; no spatial hypothesis test')
}
views<-do.call(rbind,viewmaps);keys<-do.call(rbind,viewkeys)
stopifnot(nrow(keys)==nrow(cells),!anyDuplicated(keys$cell_key),setequal(keys$cell_key,cells$cell_key))
write.table(views,file.path(sd,'spatial_view_bounds.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
gz<-gzfile(file.path(sd,'spatial_view_cell_keys.tsv.gz'),'wt');write.table(keys,gz,sep='\t',quote=FALSE,row.names=FALSE);close(gz)
write.table(do.call(rbind,dimensions),file.path(x,'review/figure_dimensions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,maps),file.path(x,'review/panel_source_map.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(x,'review/R_sessionInfo.txt'))
