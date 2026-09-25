#!/usr/bin/env Rscript
# Presentation only. Never source any original analysis script.
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(tidyr);library(cowplot);library(ggrepel)})
args <- commandArgs(TRUE)
stopifnot(length(args)>=2)
src <- normalizePath(args[1],winslash='/',mustWork=TRUE)
out <- normalizePath(args[2],winslash='/',mustWork=TRUE)
stopifnot(!startsWith(out,paste0(src,'/')),!startsWith(src,paste0(out,'/')),src!=out)
options(stringsAsFactors=FALSE,scipen=5)
for (folder in c('source_data/panel_sources','figures','supplementary','review','logs')) dir.create(file.path(out,folder),recursive=TRUE,showWarnings=FALSE)
catalog <- list(); dimensions <- list()
rd <- function(rel) {
 p<-file.path(src,rel); if(!file.exists(p))stop('Missing input: ',p)
 d<-read.delim(p,check.names=FALSE,na.strings=c('NA',''),stringsAsFactors=FALSE)
 d$.source_file<-rel; d$.source_row<-seq_len(nrow(d)); d
}
tab <- function(n) {
 p<-list.files(file.path(src,'tables_for_article'),pattern=paste0('^Table_S',n,'_'),full.names=FALSE)
 stopifnot(length(p)==1);rd(file.path('tables_for_article',p))
}
emit <- function(id,d,unit,transformation='Original estimates; presentation labels only') {
 p<-file.path(out,'source_data','panel_sources',paste0(id,'.tsv'))
 write.table(d,p,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
 sources<-if('.source_file'%in%names(d))paste(unique(d$.source_file),collapse=';')else 'See companion input rows'
 catalog[[length(catalog)+1]]<<-data.frame(panel=id,file=p,original_source=sources,record_key='source_file + source_row (1-based data row, excludes header)',unit=unit,transformation=transformation)
 invisible(d)
}
savefig <- function(p,id,h,w=7.2) {
 folder<-if(startsWith(id,'Supplementary'))'supplementary' else 'figures'
 base<-file.path(out,folder,id)
 ggsave(paste0(base,'.pdf'),p,width=w,height=h,units='in',device=cairo_pdf,bg='white')
 ggsave(paste0(base,'.tiff'),p,width=w,height=h,units='in',dpi=600,device=grDevices::tiff,compression='lzw',bg='white')
 ggsave(paste0(base,'_preview.png'),p,width=w,height=h,units='in',dpi=200,device=grDevices::png,bg='white')
 dimensions[[length(dimensions)+1]]<<-data.frame(figure=id,width_in=w,height_in=h,tiff_dpi=600)
 message('Saved ',id)
}
ink<-'#22313D';blue<-'#2C79A6';red<-'#C14A43';gold<-'#C59633';purple<-'#795AA6';grey<-'#84929A'
theme_set(theme_classic(base_size=9,base_family='sans')+theme(text=element_text(colour=ink),axis.text=element_text(colour=ink),plot.title=element_text(size=10,face='bold'),plot.subtitle=element_text(size=8),strip.background=element_rect(fill='#EDF1F3',colour=NA),strip.text=element_text(size=8,face='bold'),legend.title=element_text(size=8),legend.text=element_text(size=7.5),legend.key.size=unit(3,'mm'),plot.margin=margin(7,7,7,7)))
wrap <- function(x,n=26)vapply(x,function(s)paste(strwrap(s,n),collapse='\n'),character(1))
fmt <- function(x)ifelse(is.na(x),'NE',formatC(x,digits=3,format='g'))
labelled <- function(p,l)ggdraw(p)+draw_plot_label(l,x=0,y=1,hjust=0,vjust=1,size=12)
gridp <- function(...,labs,ncol=1,rel_heights=1,rel_widths=1)plot_grid(...,labels=labs,ncol=ncol,rel_heights=rel_heights,rel_widths=rel_widths,label_size=12,label_fontface='bold')
mods<-c('CURATED_IFN_VISIBILITY','ICI_COLITIS_EPITHELIAL_UP_STRICT','ICI_COLITIS_EPITHELIAL_UP_NON_IFN','ICI_COLITIS_EPITHELIAL_LOSS_STRICT')
mlab<-c('IFN visibility','Disease-up','Non-IFN residual','Disease-depleted')
shortmod<-function(x){m<-setNames(mlab,mods);m<-c(m,CURATED_COLON_METABOLIC_FUNCTION='Metabolic reference',CURATED_EPITHELIAL_BARRIER='Barrier reference',TNFSF12_FUNCTION_RESCUE_CORE='Reciprocity core',TNFSF12_ORGANOID_UP='TNFSF12-up',TNFSF12_ORGANOID_DOWN='TNFSF12-down');v<-unname(m[x]);v[is.na(v)]<-x[is.na(v)];v}
family_order<-c('Absorptive epithelial','Immature epithelial','Mature absorptive epithelial','Secretory epithelial','Platform samples')
family_short<-c('Absorptive','Immature','Mature absorptive','Secretory','Platform samples')
fam<-function(x)factor(x,levels=rev(family_order))
pathlab<-function(x)wrap(gsub('_',' ',sub('^HALLMARK_','',x)),26)
nesfill<-function(name='NES')scale_fill_gradient2(low=blue,mid='white',high=red,midpoint=0,na.value='#B8BDC2',limits=function(z)c(-max(abs(z)),max(abs(z))),name=name)
missing_grid<-function(d,x,y,id) {
 z<-tidyr::complete(d,!!rlang::sym(x),!!rlang::sym(y))
 z$mark[is.na(z$NES)]<-'NE'
 m<-z[is.na(z$NES),c(x,y,'NES')];m$status<-'Not estimable or not reported; not zero'
 write.table(m,file.path(out,'source_data',paste0(id,'_missing_grid.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
 z
}

# Figure 1: full-width workflow with supportive branches and explicit units.
coh<-rd('tables_for_article/Table_1_study_cohorts_and_roles.tsv')
emit('F1A_cohort_inputs',coh,'Patients/donors, cells/nuclei, released samples/spots, and ligands remain distinct')
nodes<-data.frame(x=seq(1,9,2),y=1,label=c('Receiver\nremodelling','Frozen epithelial\nprogrammes','Human organoid\nperturbation','Tumour / pathway\nguardrails','Evidence\nboundaries'),fill=c('#DDECF4','#E4EFE9','#F5ECD7','#EEE7F2','#EFE4E2'))
flow<-ggplot(nodes,aes(x,y))+geom_tile(aes(fill=fill),width=1.7,height=.7,colour='#D0D8DD')+scale_fill_identity()+geom_text(aes(label=label),size=2.8,fontface='bold',lineheight=1.1)+
 geom_segment(data=data.frame(x=seq(1.85,7.85,2),xe=seq(2.15,8.15,2)),aes(x=x,xend=xe,y=1,yend=1),inherit.aes=FALSE,arrow=arrow(length=unit(1.5,'mm')),colour=ink)+
 annotate('segment',x=1,xend=1,y=.64,yend=.20,linetype=2,colour=blue)+annotate('text',x=1,y=.04,label='TCR support',size=2.7,colour=blue)+
 annotate('segment',x=3,xend=3,y=.64,yend=.20,linetype=2,colour=grey)+annotate('text',x=3.25,y=.04,label='Lung sensitivity\n(source-confounded)',size=2.6,colour=ink)+
 coord_cartesian(xlim=c(.05,9.95),ylim=c(-.25,1.48),clip='off')+theme_void()+labs(title='From tissue programmes to perturbation selectivity')+theme(plot.title=element_text(size=11,face='bold'),plot.margin=margin(10,8,5,8))
core<-rd('02_results/core_cohort_summary.tsv');emit('F1B_core_inputs',core,'Patients and cells separately')
sp<-rd('02_results/GSE210037_spatial_pseudobulk_metadata.tsv');emit('F1B_platform_inputs',sp,'Released sample; spots are observations within a sample')
clin<-tab(24);emit('F1B_clinical_inputs',clin,'Patients by endpoint; cohorts overlap across endpoints')
land<-coh%>%mutate(row=rev(seq_len(n())),display=paste(context,cohort,sep=' | '),units=evaluable_samples)
land$units[land$cohort=='GSE228597']<-paste0(land$units[land$cohort=='GSE228597'],'; ',format(sum(core$n_cells[core$dataset=='GSE228597']),big.mark=','),' cells/nuclei')
colonsp<-sp[sp$tissue%in%c('irAE_colitis','healthy_colon'),]
land$units[land$cohort=='GSE210037']<-paste0(nrow(sp),' samples, ',format(sum(sp$n_spots),big.mark=','),' spots overall; colon: ',nrow(colonsp),' samples, ',format(sum(colonsp$n_spots),big.mark=','),' spots')
emit('F1B_display',land,'Separate patients/donors and cells/spots','Counts appended from F1B_core_inputs and F1B_platform_inputs; original evaluable_samples retained')
land$units<-wrap(land$units,64)
land$display<-wrap(land$display,25)
p1b<-ggplot(land,aes(y=row))+geom_hline(yintercept=seq(.5,nrow(land)+.5,1),colour='#E0E5E8',linewidth=.3)+geom_text(aes(x=0,label=display),hjust=0,size=2.65,fontface='bold')+geom_text(aes(x=3.1,label=units),hjust=0,size=2.6,lineheight=1.12)+coord_cartesian(xlim=c(0,10),clip='off')+theme_void()+labs(title='Cohorts and observation units')+theme(plot.title=element_text(face='bold',size=10),plot.margin=margin(12,8,10,8))
newmeta<-read.delim(file.path(out,'spatial_extension/source_data/sample_manifest.tsv'),check.names=FALSE)
newqc<-read.delim(file.path(out,'spatial_extension/source_data/section_QC.tsv'),check.names=FALSE)
newcounts<-data.frame(sections=nrow(newmeta),patients=dplyr::n_distinct(newmeta$patient),ICI_patients=dplyr::n_distinct(newmeta$patient[newmeta$condition=='ICI colitis']),healthy_patients=dplyr::n_distinct(newmeta$patient[newmeta$condition=='Healthy']),UC_patients=dplyr::n_distinct(newmeta$patient[newmeta$condition=='UC']),tissue_spots=sum(newqc$kept_tissue_spots))
write.table(newcounts,file.path(out,'source_data/F1C_new_spatial_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
xmeta<-read.delim(file.path(out,'xenium_extension/source_data/patient_compartment_manifest.tsv'),check.names=FALSE)
xmeta<-xmeta[xmeta$patient!='unassigned',]
xcounts<-data.frame(patients=dplyr::n_distinct(xmeta$patient),cells=sum(xmeta$source_cells),role='Spatial follow-through; overlapping V1 donors, not independent validation')
write.table(xcounts,file.path(out,'source_data/F1C_xenium_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
p1c<-ggdraw()+draw_label(paste0('Additional spatial cohort | GSE189184 + adult reference GSE158328\n',newcounts$sections,' sections / ',newcounts$patients,' donors: ',newcounts$ICI_patients,' ICI colitis, ',newcounts$healthy_patients,' healthy, ',newcounts$UC_patients,' UC; ',format(newcounts$tissue_spots,big.mark=','),' tissue spots\nPatient-level transfer and registered tissue maps: Supplementary Figs. 8 and 9A-F\n\nTargeted Xenium follow-through | ',xcounts$patients,' donor labels; ',format(xcounts$cells,big.mark=','),' assigned cells\nOverlapping V1 donors; cellular localization: Supplementary Figs. 10-12'),x=.025,hjust=0,size=8)
savefig(gridp(flow,p1b,p1c,labs=c('a','b','c'),rel_heights=c(1.7,5.4,1.15)),'Figure_1_study_design',8.3)

# Figure 2: original classification and all-pathway scatter; raw lung NES.
s52<-tab(52);s53<-tab(53);s54<-tab(54);s57<-tab(57)
focus<-s53%>%filter(cross_organ_class%in%c('shared_up_receiver_program','organ_divergent_receiver_program'))
paths<-focus%>%arrange(class_order,desc(heart_median_NES+colon_median_NES))%>%pull(pathway)
cts<-c('Heart | Endothelial cells','Heart | Fibroblasts','Heart | Mural cells','Colon | Endothelial','Colon | Epithelial','Colon | Mesenchymal Stromal')
heat<-s52%>%filter(pathway%in%paths)%>%mutate(compartment=factor(paste(ifelse(receiver_organ=='heart','Heart','Colon'),celltype,sep=' | '),levels=cts),path_label=factor(pathway,levels=rev(paths)),mark=ifelse(is.na(NES),'NE',ifelse(!is.na(padj)&padj<.05,'*','')))
emit('F2A',heat,'Patient-level compartment GSEA')
p2a<-ggplot(heat,aes(compartment,path_label,fill=NES))+geom_tile(colour='white')+geom_text(aes(label=mark),size=3)+nesfill()+scale_y_discrete(labels=pathlab)+scale_x_discrete(labels=c('Heart\nendothelial','Heart\nfibroblast','Heart\nmural','Colon\nendothelial','Colon\nepithelial','Colon\nstromal'))+labs(title='Shared and context-divergent programmes',subtitle='Shared and divergent classes; * compartment FDR < 0.05',x=NULL,y=NULL)+theme(axis.text.y=element_text(size=7.6),axis.text.x=element_text(size=7.3),legend.position='right')
p2a<-p2a+geom_hline(yintercept=sum(focus$cross_organ_class=='organ_divergent_receiver_program')+.5,colour=ink,linewidth=.45)
sc<-s53%>%mutate(class=case_when(cross_organ_class=='shared_up_receiver_program'~'Shared up',cross_organ_class=='organ_divergent_receiver_program'~'Context-divergent',TRUE~'Other'),label=ifelse(class=='Other','',gsub('_',' ',sub('^HALLMARK_','',pathway))))
emit('F2B',sc,'Pathway; organ medians of available compartments');emit('F2B_statistics',s54,'Original Spearman summary')
p2b<-ggplot(sc,aes(heart_median_NES,colon_median_NES,colour=class))+geom_hline(yintercept=0,colour='#C7CFD4')+geom_vline(xintercept=0,colour='#C7CFD4')+geom_point(size=1.8)+geom_text_repel(aes(label=label),size=2,max.overlaps=Inf,seed=1,box.padding=.25,show.legend=FALSE)+scale_colour_manual(values=c('Shared up'=red,'Context-divergent'=blue,'Other'=grey))+labs(title='Complete Hallmark comparison',subtitle=paste0('rho = ',sprintf('%.4f',s54$spearman_rho),'; P = ',fmt(s54$p)),x='Heart median NES',y='Colon median NES',colour=NULL)+theme(legend.position='bottom',legend.text=element_text(size=6.6))
lung<-s57%>%filter(pathway%in%mods)%>%mutate(module_label=factor(paste0(shortmod(pathway),ifelse(expected_direction_in_irAE=='down','\n(expected down)','\n(expected up)')),levels=rev(paste0(mlab,c('\n(expected up)','\n(expected up)','\n(expected up)','\n(expected down)')))),ct=gsub('Transit Amplifying','TA',gsub('Monocyte derived Macro','Monocyte-derived macro',gsub('_',' ',celltype))),sig=ifelse(is.na(padj),'NE',ifelse(padj<.05,'FDR < 0.05','FDR >= 0.05')))
emit('F2C',lung,'BALF immune-cell-type estimate; source-confounded')
p2c<-ggplot(lung,aes(NES,module_label,colour=ct,shape=sig))+geom_vline(xintercept=0,colour=grey,linetype=2)+geom_point(position=position_dodge(width=.45,orientation='y'),size=2.2)+scale_shape_manual(values=c('FDR < 0.05'=16,'FDR >= 0.05'=1,'NE'=4))+labs(title='Lung immune sensitivity',subtitle='Raw NES; source-confounded\nFilled: FDR < 0.05; open: FDR >= 0.05',x='Raw NES',y=NULL,colour=NULL,shape=NULL)+theme(legend.position='bottom',legend.text=element_text(size=6.5),axis.text.y=element_text(size=7))+guides(colour=guide_legend(ncol=1),shape=guide_legend(ncol=1))
figure2<-ggdraw()+draw_plot(p2a,0,.50,1,.50)+draw_plot(p2b,0,0,.50,.50)+draw_plot(p2c,.50,0,.50,.50)+draw_plot_label(c('a','b','c'),x=c(0,0,.50),y=c(1,.50,.50),size=12,hjust=0,vjust=1)
savefig(figure2,'Figure_2_receiver_programmes',9)

# Figure 3: freeze membership, estimates and bootstrap intervals.
s21<-tab(21);s22<-tab(22);s56<-tab(56)
val<-s22%>%filter(module%in%mods,cohort_role%in%c('independent_primary','independent_spatial_validation'))%>%mutate(oriented_NES=ifelse(expected_direction_in_irAE=='down',-NES,NES),validation=fam(ifelse(validation_dataset=='GSE210037','Platform samples',celltype)),module_label=factor(shortmod(module),levels=mlab),mark=ifelse(is.na(NES),'NE',ifelse(padj<.05&direction_concordant,'*','')))
emit('F3A',val,'Epithelial family within one external cohort; independent platform samples','Display NES = -raw NES for expected down; +raw NES otherwise')
p3a<-ggplot(val,aes(module_label,validation,fill=oriented_NES))+geom_tile(colour='white')+geom_text(aes(label=mark),size=3.5)+geom_hline(yintercept=1.5,colour=ink,linewidth=.55)+nesfill('Oriented NES')+scale_y_discrete(labels=setNames(family_short,family_order))+labs(title='Frozen programme transfer',subtitle='Positive = expected disease direction; * concordant FDR < 0.05',x=NULL,y=NULL)+theme(axis.text.x=element_text(size=8),legend.position='right')
auc<-s56%>%filter(module%in%mods,dataset%in%c('GSE206300','GSE210037'))%>%mutate(validation=fam(ifelse(dataset=='GSE210037','Platform samples',celltype)),module_label=factor(shortmod(module),levels=mlab))
emit('F3B',auc,'Patient or released sample; saved AUC and saved bootstrap 95% CI')
p3b<-ggplot(auc,aes(auc_disease_higher,validation,colour=module_label))+geom_vline(xintercept=.5,linetype=2,colour=grey)+geom_hline(yintercept=1.5,colour='#CFD7DC')+geom_segment(aes(x=bootstrap_auc_ci_low,xend=bootstrap_auc_ci_high,yend=validation),linewidth=.6)+geom_point(size=2)+facet_wrap(~module_label,ncol=2)+scale_colour_manual(values=setNames(c(red,gold,purple,blue),mlab),guide='none')+scale_x_continuous(limits=c(0,1.02),breaks=c(0,.5,1))+scale_y_discrete(labels=setNames(family_short,family_order))+labs(title='Score repeatability',subtitle='Original AUC and bootstrap intervals',x='Direction-oriented AUC',y=NULL)+theme(axis.text.y=element_text(size=7.5))
members<-s21%>%filter(module%in%mods);emit('F3C_members',members,'Gene membership; no reselection')
sizes<-members%>%distinct(module,gene)%>%count(module,name='n_genes')%>%mutate(module_label=factor(shortmod(module),levels=rev(mlab)))
emit('F3C_sizes',sizes,'Distinct genes','Count distinct module-gene pairs in F3C_members')
p3c<-ggplot(sizes,aes(n_genes,module_label,fill=module_label))+geom_col(width=.6)+geom_text(aes(label=n_genes),hjust=-.15,size=3)+scale_fill_manual(values=setNames(c(red,gold,purple,blue),mlab),guide='none')+scale_x_continuous(expand=expansion(mult=c(0,.16)))+labs(title='Unchanged module membership',x='Genes',y=NULL)
savefig(gridp(p3a,p3b,p3c,labs=c('a','b','c'),rel_heights=c(2.1,3.8,1.8)),'Figure_3_frozen_module_transfer',8)

# Figure 4: retained disease-perturbation connection, separate estimands.
s23<-tab(23);s47<-tab(47);s49<-tab(49);s50<-tab(50)
candidates<-s23%>%filter(function_rescue_gate==TRUE,!is_control)%>%pull(perturbation)
screen<-s23%>%filter(is.finite(ici_loss_median),is.finite(ici_up_non_ifn_median))%>%mutate(class=ifelse(perturbation%in%candidates,'Directional nominees','Other evaluable'),label=ifelse(perturbation%in%candidates,perturbation,''),display_suppression=-ici_up_non_ifn_median)
emit('F4A',screen,'Perturbation; median of retained cell-type NES','Y = negative saved median non-IFN NES; original candidate set unchanged')
p4a<-ggplot(screen,aes(ici_loss_median,display_suppression,colour=class))+geom_hline(yintercept=0,colour=grey)+geom_vline(xintercept=0,colour=grey)+geom_point(size=1.8)+geom_text_repel(aes(label=label),size=2.3,max.overlaps=Inf,seed=1,show.legend=FALSE)+scale_colour_manual(values=c('Directional nominees'=gold,'Other evaluable'=grey))+labs(title='ICI programme response to ligands',x='Disease-depleted programme\nmedian NES',y='Non-IFN suppression\n(-median NES)',colour=NULL)+theme(legend.position='bottom')
og<-rd('02_results/GSE313368_irAE_module_GSEA_all.tsv.gz')
tn<-og%>%filter(perturbation=='TNFSF12',module%in%mods)%>%mutate(module_label=factor(shortmod(module),levels=rev(mlab)),sig=ifelse(screen_padj<.05,'Screen FDR < 0.05','Screen FDR >= 0.05'))
emit('F4B',tn,'Cell-type NES; BH across perturbations within cell type and module')
p4b<-ggplot(tn,aes(NES,module_label,shape=sig))+geom_vline(xintercept=0,linetype=2,colour=grey)+geom_point(size=2.6,colour=blue)+scale_shape_manual(values=c('Screen FDR < 0.05'=16,'Screen FDR >= 0.05'=1))+labs(title='TNFSF12 programme effects',subtitle='Transit-amplifying-cell estimates',x='Raw NES',y=NULL,shape=NULL)+theme(legend.position='bottom',legend.text=element_text(size=7),axis.text.y=element_text(size=7.3))+guides(shape=guide_legend(ncol=1))
ef<-rd('02_results/GSE313368_celltype_module_effects.tsv.gz')%>%filter(perturbation=='TNFSF12',analysis_scope=='perturbation_vs_cytomix',module%in%c('IBD1','IBD2','GP20'))%>%mutate(sig=ifelse(padj<.05,'Publisher FDR < 0.05','Publisher FDR >= 0.05'),ct=gsub('Transit Amplifying','TA',gsub('Monocyte derived Macro','Monocyte-derived macro',gsub('_',' ',celltype))))
emit('F4C',ef,'Published cell-type Cohen d; donor/well aggregation performed by source authors')
p4c<-ggplot(ef,aes(effect_size,module,colour=ct,shape=sig))+geom_vline(xintercept=0,colour=grey,linetype=2)+geom_point(position=position_dodge(width=.5,orientation='y'),size=2.3)+scale_shape_manual(values=c('Publisher FDR < 0.05'=16,'Publisher FDR >= 0.05'=1))+labs(title='Published programme context',subtitle='TNFSF12 versus matched cytomix control',x="Published Cohen's d",y=NULL,colour=NULL,shape=NULL)+theme(legend.position='bottom')+guides(colour=guide_legend(ncol=1),shape=guide_legend(ncol=1))
off<-s47%>%filter(perturbation=='TNFSF12')%>%mutate(path_label=factor(pathway,levels=rev(unique(pathway))),mark=ifelse(is.na(NES),'NE',paste0(sprintf('%.2f',NES),ifelse(!is.na(padj)&padj<.05,'*',''))))
emit('F4D',off,'Cell-type pathway estimate; original Hallmark-family FDR')
p4d<-ggplot(off,aes(celltype,path_label,fill=NES))+geom_tile(colour='white')+geom_text(aes(label=mark),size=2.5)+nesfill()+scale_y_discrete(labels=pathlab)+scale_x_discrete(labels=function(x)gsub('_',' ',x))+labs(title='Broader pathway response',subtitle='Raw NES; * original Hallmark FDR < 0.05',x=NULL,y=NULL)+theme(axis.text.y=element_text(size=6.8),axis.text.x=element_text(size=7.2),legend.position='bottom')
p4d<-p4d+labs(subtitle='Raw NES; * Hallmark FDR < 0.05')
rec<-s49%>%mutate(context=ifelse(cohort_role=='derivation_cohort_reciprocity_check','Discovery self-check',ifelse(dataset=='GSE210037','Platform samples',celltype)),context=factor(context,levels=c('Discovery self-check',family_order)),signature=factor(shortmod(pathway),levels=c('TNFSF12-up','TNFSF12-down','Reciprocity core')),mark=ifelse(is.na(NES),'NE',paste0(sprintf('%.2f',NES),ifelse(!is.na(padj)&padj<.05,'*',''))))
emit('F4E',rec,'Epithelial family / sample-level contrast; discovery is not independent validation')
emit('F4E_summary',s50,'Saved cohort-role aggregate with original denominators')
p4e<-ggplot(rec,aes(context,signature,fill=NES))+geom_tile(colour='white')+geom_text(aes(label=mark),size=2.7)+geom_vline(xintercept=c(1.5,5.5),colour=ink,linewidth=.5)+nesfill()+scale_x_discrete(labels=c('Discovery\nself-check','Absorptive','Immature','Mature\nabsorptive','Secretory','Platform\nsamples'))+labs(title='Epithelial reciprocity across disease contexts',subtitle='Raw disease NES; * FDR < 0.05; core fixed before external transfer',x=NULL,y=NULL)+theme(axis.text.x=element_text(size=7.2),axis.text.y=element_text(size=8),legend.position='right')
# Explicitly incorporate the completed, separately executed donor analysis.
# The existing V1 b/c panels remain available as Supplementary Figure 14.
orun<-jsonlite::fromJSON(file.path(out,'organoid_donor_extension/review/analysis_complete.json'))
stopifnot(orun$status=='PASS')
orgread<-function(n){p<-file.path(orun$output,n);z<-read.delim(p,check.names=FALSE);z$.source_file<-p;z$.source_row<-seq_len(nrow(z));z}
ng<-orgread('module_GSEA.tsv');ns<-orgread('donor_module_scores.tsv')
old4b<-p4b+labs(title='V1 export-based TNFSF12 estimates')
old4c<-p4c
ref<-ng%>%filter(ligand=='No_Cytomix')%>%mutate(module_label=factor(shortmod(module),levels=rev(mlab)),ct=factor(celltype,levels=c('Stem','Transit_Amplifying','Colonocyte','Goblet')),mark=ifelse(is.na(NES),'NE',paste0(sprintf('%.2f',NES),ifelse(FDR_predeclared_family<.05,'*',''))))
emit('SF14C_REFERENCE',ref,'Gene-set NES from donor-adjusted No_Cytomix vs cytomix model; separate 16-test BH family','No new candidate nomination')
pref<-ggplot(ref,aes(ct,module_label,fill=NES))+geom_tile(colour='white')+geom_text(aes(label=mark),size=3)+nesfill()+scale_x_discrete(labels=c('Stem','TA','Colonocyte','Goblet'))+labs(title='Assay reference without cytomix',subtitle='Raw NES; * separate 16-test FDR < 0.05',x=NULL,y=NULL)
savefig(gridp(gridp(old4b,old4c,labs=c('a','b'),ncol=2),pref,labs=c('','c'),rel_heights=c(4,2.3)),'Supplementary_Figure_14_original_and_assay_context',6.8)
n4<-ng%>%filter(ligand=='TNFSF12')%>%mutate(module_label=factor(shortmod(module),levels=rev(mlab)),ct=factor(celltype,levels=c('Stem','Transit_Amplifying','Colonocyte','Goblet')),sig=ifelse(FDR_predeclared_family<.05,'FDR < 0.05','FDR >= 0.05'))
emit('F4B_DONOR',n4,'Gene-set NES from three-donor paired DESeq2 Wald ranks; BH across 16 primary tests','Separate new analysis, not a rerun of publisher effect ratios')
p4a<-p4a+labs(title='V1 directional nomination')
p4b<-ggplot(n4,aes(NES,module_label,colour=ct,shape=sig))+geom_vline(xintercept=0,linetype=2,colour=grey)+geom_point(position=position_dodge(width=.6,orientation='y'),size=2.2)+scale_colour_manual(values=c(blue,red,purple,gold),labels=c('Stem','TA','Colonocyte','Goblet'))+scale_shape_manual(values=c('FDR < 0.05'=16,'FDR >= 0.05'=1))+labs(title='Donor-adjusted TNFSF12 response',subtitle='Fixed ICI sets; three paired donors',x='Wald-rank GSEA NES',y=NULL,colour=NULL,shape=NULL)+theme(legend.position='bottom')+guides(colour=guide_legend(ncol=2),shape=guide_legend(ncol=1))
delta<-ns%>%filter(ligand=='TNFSF12',module=='ICI_COLITIS_EPITHELIAL_LOSS_STRICT')%>%select(celltype,donor,treatment,raw_mean_rank)%>%pivot_wider(names_from=treatment,values_from=raw_mean_rank)%>%mutate(change=Treatment-Control,ct=factor(celltype,levels=rev(c('Stem','Transit_Amplifying','Colonocyte','Goblet'))))
delta$.source_file<-file.path(orun$output,'donor_module_scores.tsv');delta$.source_row<-NA_integer_
emit('F4C_DONOR',delta,'Donor; paired programme-rank change','Keys ligand=TNFSF12 + module + celltype + donor + treatment; Treatment minus Control')
p4c<-ggplot(delta,aes(change,ct,colour=donor))+geom_vline(xintercept=0,linetype=2,colour=grey)+geom_point(position=position_dodge(width=.5,orientation='y'),size=2.3)+scale_y_discrete(labels=c('Stem'='Stem','Transit_Amplifying'='TA','Colonocyte'='Colonocyte','Goblet'='Goblet'))+labs(title='Disease-depleted programme by donor',subtitle='Each point is one paired donor\nPositive = increased relative expression',x='Treatment minus control mean rank',y=NULL,colour='Donor')+theme(legend.position='bottom')
p4d<-p4d+labs(subtitle='V1 export-based NES\n* original FDR < 0.05')
p4e<-p4e+labs(subtitle='Original disease NES; * FDR < 0.05; V1 gene sets remain frozen')
savefig(gridp(gridp(p4a,p4b,labs=c('a','b'),ncol=2),gridp(p4c,p4d,labs=c('c','d'),ncol=2,rel_widths=c(.95,1.05)),p4e,labs=c('','','e'),rel_heights=c(3,4.5,2)),'Figure_4_organoid_reciprocity',9.8)

# Figure 5: all feature categories, with no outcome-dependent feature choice.
s27<-tab(27);s29<-tab(29);s43<-tab(43);s25<-tab(25)
category<-function(x)case_when(x%in%c('direct_ligand_gene','direct_receptor_gene')~'GENE',x=='ligand_receptor_axis'~'AXIS',TRUE~'SIGNATURE')
feature_label<-function(x)gsub('_COLONOCYTE',' COL',gsub('_TA',' TA',gsub('TNFSF12_ICI_LOSS_RESCUE_CORE','TNFSF12 core',gsub('^GENE_|^AXIS_|^PERTSIG_','',x))))
resp<-s27%>%mutate(category=factor(category(feature_type),levels=c('GENE','AXIS','SIGNATURE')),feature=feature_label(feature_id),feature=factor(feature,levels=rev(unique(feature))))
emit('F5A',resp,'Original response meta-estimate and 95% CI; all feature classes')
p5a<-ggplot(resp,aes(estimate,feature))+geom_vline(xintercept=0,linetype=2,colour=grey)+geom_segment(aes(x=ci_low,xend=ci_high,yend=feature),colour=blue,linewidth=.45)+geom_point(aes(shape=fdr<.05),size=1.8,colour=blue)+facet_grid(category~.,scales='free_y',space='free_y')+scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),labels=c('FDR >= 0.05','FDR < 0.05'))+labs(title='Tumour response',subtitle="All original feature categories",x="Hedges' g (higher in responders > 0)",y=NULL,shape=NULL)+theme(axis.text.y=element_text(size=6.8),legend.position='bottom',strip.text.y=element_text(angle=90,size=7))
surv<-s29%>%filter(meta_scope=='OS_only')%>%mutate(category=factor(category(feature_type),levels=c('GENE','AXIS','SIGNATURE')),feature=feature_label(feature_id),feature=factor(feature,levels=levels(resp$feature)))
if(!nrow(surv))stop('No OS_only estimates')
emit('F5B',surv,'Original OS hazard ratio and 95% CI; no Cox or pooling rerun')
p5b<-ggplot(surv,aes(hazard_ratio,feature))+geom_vline(xintercept=1,linetype=2,colour=grey)+geom_segment(aes(x=ci_low_hr,xend=ci_high_hr,yend=feature),colour=purple,linewidth=.45)+geom_point(aes(shape=fdr<.05),size=1.8,colour=purple)+facet_grid(category~.,scales='free_y',space='free_y')+scale_x_log10(breaks=c(.5,.75,1,1.5,2))+ scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),labels=c('FDR >= 0.05','FDR < 0.05'))+labs(title='Overall survival',subtitle='Pooled associations per expression SD',x='Hazard ratio (log scale)',y=NULL,shape=NULL)+theme(axis.text.y=element_text(size=6.8),legend.position='bottom',strip.text.y=element_text(angle=90,size=7))
ind<-rd('02_results/independent_colon_patient_level_DE_all.tsv.gz')%>%filter(feature_level=='family',celltype%in%family_order,contrast=='irColitis_vs_ICI_control_PD1_only',sub('^.*\\|','',gene)=='SORL1')
emit('F5C_family',ind,'SORL1 original family-specific DE')
sorl<-s43%>%filter(gene=='SORL1');emit('F5C_summary',sorl,'Saved SORL1 results; no minimum-family FDR attached to median')
new_sor<-read.delim(gzfile(file.path(orun$output,'gene_results/TNFSF12__Transit_Amplifying.tsv.gz')))%>%mutate(.source_file=file.path(orun$output,'gene_results/TNFSF12__Transit_Amplifying.tsv.gz'),.source_row=row_number())%>%filter(gene=='SORL1');stopifnot(nrow(new_sor)==1)
emit('F5C_DONOR',new_sor,'SORL1 gene-specific paired donor estimate; FDR over genes in its contrast','Added alongside unchanged V1 export estimate; not AXIS or SIGNATURE significance')
st<-bind_rows(data.frame(context='Organoid TA',estimate=sorl$organoid_log2FC,q=sorl$organoid_fdr),data.frame(context='Organoid TA donor model',estimate=new_sor$log2FoldChange,q=new_sor$padj),data.frame(context=as.character(ind$celltype),estimate=ind$logFC,q=ind$FDR),data.frame(context='Platform samples',estimate=sorl$spatial_logFC,q=sorl$spatial_fdr))%>%mutate(context=factor(context,levels=rev(c('Organoid TA','Organoid TA donor model',family_order))),label=paste0('q=',fmt(q)))
emit('F5C_display',st,'Cell-type / platform log2 fold change','Values from F5C_family and F5C_summary; original direction, no fabricated CI')
p5c<-ggplot(st,aes(estimate,context))+geom_vline(xintercept=0,colour=grey,linetype=2)+geom_point(colour=purple,size=2.3)+geom_text(aes(label=label),hjust=-.14,size=2.6)+scale_y_discrete(labels=c('Platform samples'='Platform samples','Secretory epithelial'='Secretory','Mature absorptive epithelial'='Mature absorptive','Absorptive epithelial'='Absorptive','Immature epithelial'='Immature','Organoid TA'='Organoid TA / V1','Organoid TA donor model'='Organoid TA / paired'))+scale_x_continuous(expand=expansion(mult=c(.08,.35)))+labs(title='SORL1 across epithelial contexts',subtitle='Gene-specific estimates and FDR',x='log2 fold change',y=NULL)
coverage<-clin%>%select(cohort,n_response_evaluable,n_survival_evaluable,survival_endpoint,.source_file,.source_row);emit('F5D',coverage,'Endpoint-specific patients; no risk imputation for missing endpoint')
ct<-clin%>%mutate(label=paste0(cohort,'\nResponse n=',n_response_evaluable,'; ',ifelse(is.na(survival_endpoint),'no survival',paste0(survival_endpoint,' n=',n_survival_evaluable))),row=rev(seq_len(n())))
p5d<-ggplot(ct,aes(0,row,label=label))+geom_text(hjust=0,size=2.8,lineheight=1.1)+coord_cartesian(xlim=c(0,1),ylim=c(-.6,6.1),clip='off')+theme_void()+labs(title='Endpoint coverage')+theme(plot.title=element_text(size=10,face='bold'),plot.margin=margin(14,7,7,19))+annotate('text',x=0,y=-.2,hjust=0,label='Associations do not establish\nintervention benefit, harm or safety.',size=2.7)
savefig(gridp(gridp(p5a,p5b,labs=c('a','b'),ncol=2),gridp(p5c,p5d,labs=c('c','d'),ncol=2,rel_widths=c(1.25,.75)),labs=c('',''),rel_heights=c(6.7,2.8)),'Figure_5_tumour_guardrails',9.8)

# Figure 6: patient points and saved statistics only; no fitted line or band.
s31<-tab(31);s35<-tab(35);s36<-tab(36);s37<-tab(37)
metrics<-c('expanded_cell_fraction_ge2','normalized_clonality','cross_cd4_cd8_cell_fraction')
mt<-c('Expanded-cell\nfraction','Normalized\nclonality','CD4-CD8 shared\nclone-cell fraction');names(mt)<-metrics
tcr<-s31%>%filter(disease%in%c('CPI_colitis','HC'))%>%pivot_longer(all_of(metrics),names_to='metric',values_to='value')%>%mutate(metric_label=factor(mt[metric],levels=mt),group=ifelse(disease=='HC','Healthy','ICI colitis'))
emit('F6A',tcr,'Patient; fractions and normalized clonality are dimensionless')
ts<-s35%>%filter(analysis=='colon_CPI_colitis_vs_HC_patient_level',metric%in%metrics)%>%mutate(metric_label=factor(mt[metric],levels=mt),label=paste0('P=',fmt(p),'\nFDR=',fmt(fdr_within_analysis)))
emit('F6A_stats',ts,'Original patient-level Mann-Whitney tests')
p6a<-ggplot(tcr,aes(group,value,fill=group))+geom_boxplot(width=.48,outlier.shape=NA,alpha=.35)+geom_point(position=position_jitter(width=.07,height=0,seed=1),shape=21,size=2)+facet_wrap(~metric_label,nrow=1,scales='free_y')+geom_text(data=ts,aes(x=1.5,y=Inf,label=label),inherit.aes=FALSE,vjust=1.1,size=2.5)+scale_y_continuous(expand=expansion(mult=c(.08,.28)))+scale_fill_manual(values=c('Healthy'='#C4CCD0','ICI colitis'=red),guide='none')+labs(title='Clonal expansion in ICI colitis',subtitle=paste0(n_distinct(tcr$patient_id[tcr$group=='ICI colitis']),' ICI-colitis and ',n_distinct(tcr$patient_id[tcr$group=='Healthy']),' healthy patients; independent y-axis ranges'),x=NULL,y='Fraction / normalized index')+theme(axis.text.x=element_text(angle=20,hjust=1))
couple<-s36%>%filter(disease=='CPI_colitis');emit('F6B',couple,'Patient; raw mean within-sample gene-expression rank')
corr<-s37%>%filter(tcr_metric=='cross_cd4_cd8_cell_fraction',epithelial_module=='ICI_COLITIS_EPITHELIAL_LOSS_STRICT');emit('F6B_stats',corr,'Original Spearman correlation, P and FDR')
p6b<-ggplot(couple,aes(cross_cd4_cd8_cell_fraction,ICI_COLITIS_EPITHELIAL_LOSS_STRICT))+geom_point(colour=blue,size=2.4)+geom_text_repel(aes(label=patient_id),size=2.5,seed=1,max.overlaps=Inf)+labs(title='Exploratory epithelial association',subtitle=paste0('rho=',fmt(corr$spearman_rho),'; P=',fmt(corr$p),'; FDR=',fmt(corr$fdr)),x='CD4-CD8 shared-clone cell fraction',y='Disease-depleted programme\nraw mean expression rank')
rare<-s35%>%filter(analysis=='colon_CPI_colitis_vs_HC_patient_level',metric%in%c('rarefied_clonality_200','rarefied_top10_fraction_200'))%>%mutate(label=ifelse(metric=='rarefied_clonality_200','Rarefied clonality','Rarefied top-10 fraction'),note=paste0('P=',fmt(p),'; FDR=',fmt(fdr_within_analysis)))
emit('F6C',rare,'Original median patient difference after saved rarefaction; no new rarefaction')
p6c<-ggplot(rare,aes(median_difference,label))+geom_vline(xintercept=0,linetype=2,colour=grey)+geom_point(size=2.5,colour=gold)+geom_text(aes(label=note),hjust=-.1,vjust=-.9,size=2.7)+scale_x_continuous(expand=expansion(mult=c(.15,1.25)))+labs(title='Depth-controlled sensitivity',x='Median difference: ICI colitis - healthy',y=NULL)
savefig(gridp(p6a,p6b,p6c,labs=c('a','b','c'),rel_heights=c(3.2,3,1.7)),'Figure_6_TCR_support',8.2)

# Six supplementary figures retain the original evidence subjects.
full<-s52%>%mutate(compartment=factor(paste(ifelse(receiver_organ=='heart','Heart','Colon'),celltype,sep=' | '),levels=cts),path_label=factor(pathway,levels=rev(s53$pathway)),mark=ifelse(is.na(NES),'NE',ifelse(padj<.05,'*','')));emit('SF1',full,'Original Hallmark estimates for all receiver compartments')
full<-missing_grid(full,'compartment','path_label','SF1')
sf1<-ggplot(full,aes(compartment,path_label,fill=NES))+geom_tile(colour='white')+geom_text(aes(label=mark),size=2)+nesfill()+scale_y_discrete(labels=function(x)gsub('_',' ',sub('HALLMARK_','',x)))+scale_x_discrete(labels=c('Heart\nendothelial','Heart\nfibroblast','Heart\nmural','Colon\nendothelial','Colon\nepithelial','Colon\nstromal'))+labs(title='Complete receiver Hallmark landscape',subtitle='* FDR < 0.05; original descriptive order',x=NULL,y=NULL)+theme(axis.text.y=element_text(size=6.5),axis.text.x=element_text(size=7))
savefig(sf1,'Supplementary_Figure_1_complete_Hallmarks',10.2)
av<-s22%>%mutate(oriented_NES=ifelse(expected_direction_in_irAE=='down',-NES,NES),label=paste(validation_dataset,celltype,contrast,sep=' | '),module_label=shortmod(module),mark=ifelse(is.na(NES),'NE',ifelse(padj<.05&direction_concordant,'*','')));emit('SF2',av,'All saved contrasts; discovery and sensitivity labels explicit','Orientation only; saved values retained')
av$label<-wrap(av$label,64)
av<-missing_grid(av,'module_label','label','SF2')
sf2<-ggplot(av,aes(module_label,label,fill=oriented_NES))+geom_tile(colour='white')+geom_text(aes(label=mark),size=2.5)+nesfill('Oriented NES')+labs(title='Complete frozen-module transfer',subtitle='* concordant FDR < 0.05; discovery rows are self-checks',x=NULL,y=NULL)+theme(axis.text.y=element_text(size=6),axis.text.x=element_text(size=7,angle=35,hjust=1),legend.position='bottom')
savefig(sf2,'Supplementary_Figure_2_all_module_transfer',10)
lungall<-s57%>%mutate(module_label=paste0(shortmod(pathway),' [',expected_direction_in_irAE,']'),ct=gsub('Transit Amplifying','TA',gsub('Monocyte derived Macro','Monocyte-derived macro',gsub('_',' ',celltype))),mark=paste0(sprintf('%.2f',NES),ifelse(padj<.05,'*','')));emit('SF3',lungall,'BALF immune cell-type estimate; bracket = expected disease direction')
sf3<-ggplot(lungall,aes(ct,module_label,fill=NES))+geom_tile(colour='white')+geom_text(aes(label=mark),size=3)+nesfill()+labs(title='Source-confounded lung sensitivity',subtitle='Raw NES; * FDR < 0.05; brackets give expected disease direction',x=NULL,y=NULL)+theme(axis.text.x=element_text(angle=25,hjust=1,size=8))
savefig(sf3,'Supplementary_Figure_3_lung_sensitivity',4.8)
alloff<-s47%>%filter(perturbation%in%candidates)%>%mutate(column=paste(perturbation,gsub('Transit_Amplifying','TA',celltype),sep=' / '),path_label=pathlab(pathway),mark=ifelse(is.na(NES),'NE',ifelse(padj<.05,'*','')));emit('SF4',alloff,'Original candidate-by-cell-type off-target NES; all tested rows, not significant-only maximum')
alloff<-missing_grid(alloff,'column','path_label','SF4')
sf4<-ggplot(alloff,aes(column,path_label,fill=NES))+geom_tile(colour='white')+geom_text(aes(label=mark),size=2.5)+nesfill()+labs(title='Complete candidate off-target estimates',subtitle='* original FDR < 0.05; empty positions were not evaluated',x=NULL,y=NULL)+theme(axis.text.x=element_text(angle=55,hjust=1,size=7),axis.text.y=element_text(size=7),legend.position='bottom')
savefig(sf4,'Supplementary_Figure_4_offtarget_contexts',7.5)
s34<-tab(34)%>%filter(clone_definition=='nucleotide_exact');emit('SF5',s34,'Matched donor; exact nucleotide repertoire overlap')
sf5<-ggplot(s34,aes(compartment_b,morisita_horn_abundance,group=donor,colour=donor))+geom_line(alpha=.6)+geom_point(size=2.5)+labs(title='Paired heart-repertoire overlap',subtitle=paste0(n_distinct(s34$donor),' donors; descriptive matched comparison'),x='Comparator compartment',y='Morisita-Horn abundance overlap',colour='Donor')
savefig(sf5,'Supplementary_Figure_5_paired_repertoires',4.8)
ranking<-s43%>%filter(is.finite(independent_median_logFC),is.finite(spatial_logFC))%>%mutate(depletion=-independent_median_logFC,platform_depletion=-spatial_logFC,highlight=evidence_tier%in%c('Tier_A_convergent_downstream_candidate','Tier_B_supportive_downstream_candidate'),label=ifelse(highlight,gene,''));emit('SF6',ranking,'Gene; saved ranking, not independent target validation','Axes negate saved tissue logFC; no reranking')
sf6<-ggplot(ranking,aes(depletion,platform_depletion,colour=highlight,size=pmax(organoid_log2FC,0)))+geom_vline(xintercept=0,colour=grey)+geom_hline(yintercept=0,colour=grey)+geom_point(alpha=.75)+geom_text_repel(aes(label=label),size=2.6,seed=1,max.overlaps=Inf,show.legend=FALSE)+scale_colour_manual(values=c('TRUE'=purple,'FALSE'='#A6B0B6'),guide='none')+scale_size_continuous(range=c(1.5,4),name='Organoid log2FC')+labs(title='Downstream reciprocity-core ranking',subtitle='Original selection and estimates; exploratory markers',x='Independent epithelial depletion (-median logFC)',y='Platform depletion (-logFC)')
savefig(sf6,'Supplementary_Figure_6_downstream_ranking',5.3)

# Complete gene-level evidence, without selecting genes by their external P value.
coredefs<-tab(38)
coregenes<-sort(unique(s43$gene))
allind<-rd('02_results/independent_colon_patient_level_DE_all.tsv.gz')%>%filter(feature_level=='family',celltype%in%family_order,contrast=='irColitis_vs_ICI_control_PD1_only')%>%mutate(symbol=sub('^.*\\|','',gene))%>%filter(symbol%in%coregenes)
emit('SF7_family_inputs',allind,'Existing family-level patient pseudobulk DE; all saved core genes')
emit('SF7_summary_inputs',s43,'Existing organoid and platform gene effects; no ranking or significance selection')
geneeffects<-bind_rows(
 s43%>%transmute(gene,context='Organoid TA',effect=organoid_log2FC,q=organoid_fdr,.source_file,.source_row),
 allind%>%transmute(gene=symbol,context=celltype,effect=logFC,q=FDR,.source_file,.source_row),
 s43%>%transmute(gene,context='Platform samples',effect=spatial_logFC,q=spatial_fdr,.source_file,.source_row))
stopifnot(!anyDuplicated(geneeffects[c('gene','context')]))
emit('SF7_effects',geneeffects,'Gene-specific log2 fold change; distinct source FDR families','Organoid: ligand vs cytomix; tissues: disease vs control. Columns are not donor replicates.')
geneeffects<-tidyr::complete(geneeffects,gene=coregenes,context=c('Organoid TA',family_order))%>%mutate(context=factor(context,levels=c('Organoid TA',family_order)),gene=factor(gene,levels=rev(coregenes)),mark=ifelse(is.na(effect),'NE',ifelse(!is.na(q)&q<.05,'*','')),half=ifelse(as.character(gene)%in%head(coregenes,ceiling(length(coregenes)/2)),'a','b'))
make_core<-function(h)ggplot(geneeffects%>%filter(half==h),aes(context,gene,fill=effect))+geom_tile(colour='white',linewidth=.25)+geom_text(aes(label=mark),size=2.4)+scale_fill_gradient2(low=blue,mid='white',high=red,midpoint=0,na.value='#B8BDC2',limits=range(geneeffects$effect,na.rm=TRUE),name='log2FC')+scale_x_discrete(labels=c('Organoid TA','Absorptive','Immature','Mature absorptive','Secretory','Platform'))+labs(x=NULL,y=NULL)+theme(axis.text.y=element_text(size=7.5),axis.text.x=element_text(angle=45,hjust=1,size=7.3),legend.position='bottom',plot.margin=margin(7,12,7,7))
sf7<-gridp(make_core('a'),make_core('b'),labs=c('a','b'),ncol=2)
sf7<-plot_grid(ggdraw()+draw_label('Gene-level effects across the complete reciprocity core',fontface='bold',size=11),sf7,ggdraw()+draw_label('All members; alphabetical order. * source-specific FDR < 0.05; NE = missing.\nOrganoid: ligand vs cytomix. Tissue: disease vs control. No new model or pooled gene test.',size=8),ncol=1,rel_heights=c(.45,8,.7))
savefig(sf7,'Supplementary_Figure_7_complete_core_effects',9.2)

# Horizontal grouping must never displace the numerical x coordinate.
horizontal_checks<-lapply(c('p2c','p4b','p4c'),function(nm){
 p<-get(nm);b<-ggplot_build(p);point_layers<-which(vapply(p$layers,function(z)inherits(z$geom,'GeomPoint'),logical(1)))
 stopifnot(length(point_layers)>0)
 for(i in point_layers){
  xdata<-p$layers[[i]]$data;if(inherits(xdata,'waiver'))xdata<-p$data
  mapping<-p$mapping;if(!is.null(p$layers[[i]]$mapping$x))mapping$x<-p$layers[[i]]$mapping$x
  expected<-rlang::eval_tidy(mapping$x,data=xdata)
  actual<-b$data[[i]]$x
  stopifnot(length(expected)==length(actual),max(abs(sort(expected)-sort(actual)),na.rm=TRUE)<1e-12)
 }
 data.frame(plot=nm,point_layers=length(point_layers),numerical_x='PASS; unchanged by grouping')
})
write.table(bind_rows(horizontal_checks),file.path(out,'review','horizontal_coordinate_verification.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(bind_rows(catalog),file.path(out,'review','source_map.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(bind_rows(dimensions),file.path(out,'review','figure_dimensions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'logs','presentation_R_sessionInfo.txt'))
message('Plot files written: six main and seven supplementary figures; runtime success requires a zero process exit code.')

