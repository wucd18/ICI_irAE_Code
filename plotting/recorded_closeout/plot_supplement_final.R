#!/usr/bin/env Rscript
# Presentation only: every input is an existing table; no model/package installation.
options(warn=1);args<-commandArgs(TRUE);stopifnot(length(args)==1);O<-normalizePath(args[1],winslash='/');I<-file.path(O,'input_snapshot')
source(file.path(O,'presentation_code/csv_io.R'))
source(file.path(O,'presentation_code/theme_final.R'))
rd<-function(p)read.csv(file.path(I,p),check.names=FALSE,na.strings=c('NA',''),stringsAsFactors=FALSE)
fp<-function(p)rd(paste0('source_data/final_panels/',p,'.csv'))
old<-function(p)rd(paste0('source_data/panel_sources/',p,'.csv'))
tab<-function(n)rd(paste0('supplementary/tables/',list.files(file.path(I,'supplementary/tables'),pattern=paste0('^Table_S',n,'_.*\\.csv$'))))
for(x in c('figures','review','supplementary'))dir.create(file.path(O,x),showWarnings=FALSE)
dims<-list();plots<-list()
savep<-function(p,id,h){
 ids<-c(Supplementary_Figure_1_complete_Hallmarks='S1',Supplementary_Figure_2_all_module_transfer='S2A',Supplementary_Figure_3_lung_sensitivity='S2B',Supplementary_Figure_5_paired_repertoires='S3',Supplementary_Figure_7_complete_core_effects='S4A',Supplementary_Figure_14_original_and_assay_context='S5C')
 if(!id%in%names(ids))return(invisible(NULL))
 if(id=='Supplementary_Figure_14_original_and_assay_context'){p<-p3;h<-115}
 b<-file.path(O,'supplement/plates',ids[[id]])
 ggsave(paste0(b,'.pdf'),p,width=180,height=h,units='mm',device=cairo_pdf,bg='white')
 ggsave(paste0(b,'.png'),p,width=180,height=h,units='mm',dpi=300,bg='white')
 message('Supplement plate ',ids[[id]])
}
flabs<-c('Absorptive','Immature','Mature absorptive','Secretory','GSE210037 sample pseudobulk')
theme_update(axis.text=element_text(size=7),legend.text=element_text(size=7))
# Existing results only. Full tables and original selection retained.
raw<-tab(52);summ<-tab(53);raw$compartment<-paste(raw$receiver_organ,raw$celltype,sep=' / ');raw$pathway<-factor(raw$pathway,levels=rev(summ$pathway));raw$mark<-ifelse(is.na(raw$NES),'NE',ifelse(raw$padj<.05,'*',''))
p<-ggplot(raw,aes(compartment,pathway,fill=NES))+geom_tile(colour='white')+geom_text(aes(label=mark),size=2.3,na.rm=TRUE)+nes_scale(raw$NES)+scale_y_discrete(labels=pathlab)+scale_x_discrete(labels=function(x)wrap(gsub('_',' ',x),15))+labs(title='Complete receiver Hallmark landscape',subtitle='Original compartment estimates; * FDR < 0.05',x=NULL,y=NULL)+theme(axis.text.y=element_text(size=7))
savep(p,'Supplementary_Figure_1_complete_Hallmarks',220)
d<-old('SF2');d$module_label<-gsub('Non-IFN residual','IFN-gene-excluded residual',d$module_label);d$mark[is.na(d$mark)&!is.na(d$NES)]<-'';d$label<-gsub('Platform samples','GSE210037 sample pseudobulk',d$label,fixed=TRUE);d<-tidyr::complete(d,module_label,label);d$mark[is.na(d$NES)]<-'NE'
p<-ggplot(d,aes(module_label,label,fill=oriented_NES))+geom_tile(colour='white')+geom_text(aes(label=mark),size=2.3,na.rm=TRUE)+nes_scale(d$oriented_NES,'Oriented NES')+scale_y_discrete(labels=function(x)wrap(x,58))+labs(title='Complete fixed-programme transfer',subtitle='Positive = expected disease direction; discovery rows are self-checks',x=NULL,y=NULL)+theme(axis.text.y=element_text(size=7),axis.text.x=element_text(angle=35,hjust=1),legend.position='bottom')
savep(p,'Supplementary_Figure_2_all_module_transfer',185)
d<-tab(57);d$label<-paste0(ifelse(d$pathway%in%modids,modlab(d$pathway),d$pathway),' [',d$expected_direction_in_irAE,']');d$label<-gsub('CURATED_EPITHELIAL_BARRIER','Barrier reference',d$label);d$label<-gsub('CURATED_COLON_METABOLIC_FUNCTION','Metabolic reference',d$label);d$mark<-ifelse(is.na(d$NES),'NE',paste0(sprintf('%.2f',d$NES),ifelse(d$padj<.05,'*','')))
p<-ggplot(d,aes(celltype,label,fill=NES))+geom_tile(colour='white')+geom_text(aes(label=mark),size=2.6)+nes_scale(d$NES)+scale_x_discrete(labels=function(x)wrap(gsub('_',' ',x),15))+labs(title='Source-confounded lung module sensitivity',subtitle='Raw NES; brackets give expected disease direction; * original FDR < 0.05',x=NULL,y=NULL)
savep(p,'Supplementary_Figure_3_lung_sensitivity',125)
d<-tab(47);d<-d[d$perturbation!='NoCytomix',];d$column<-paste(d$perturbation,gsub('Transit_Amplifying','TA',d$celltype),sep=' / ');d$mark<-ifelse(is.na(d$NES),'NE',ifelse(d$padj<.05,'*',''));d<-tidyr::complete(d,column,pathway);d$mark[is.na(d$NES)]<-'NE'
p<-ggplot(d,aes(column,pathway,fill=NES))+geom_tile(colour='white')+geom_text(aes(label=mark),size=2.4,na.rm=TRUE)+nes_scale(d$NES)+scale_y_discrete(labels=pathlab)+labs(title='Published-export comparison of prior ligand cases',subtitle='Earlier export-rank Hallmark results; * source-specific FDR < 0.05; NE = unavailable',x=NULL,y=NULL)+theme(axis.text.x=element_text(angle=55,hjust=1),legend.position='bottom')
savep(p,'Supplementary_Figure_4_offtarget_contexts',185)
d<-old('SF5');p<-ggplot(d,aes(compartment_b,morisita_horn_abundance,colour=donor,group=donor))+geom_line(linewidth=.3)+geom_point(size=2)+labs(title='Matched heart-repertoire overlap',subtitle=paste0(length(unique(d$donor)),' donors; descriptive exact nucleotide-clone comparison'),x='Comparator compartment',y='Morisita-Horn abundance overlap',colour='Donor')
savep(p,'Supplementary_Figure_5_paired_repertoires',120)
d<-old('SF6');p<-ggplot(d,aes(depletion,platform_depletion,size=pmax(organoid_log2FC,0),colour=highlight))+geom_vline(xintercept=0,colour=grey,linewidth=.25)+geom_hline(yintercept=0,colour=grey,linewidth=.25)+geom_point()+geom_text_repel(aes(label=label),size=2.5,seed=1,max.overlaps=Inf,show.legend=FALSE,na.rm=TRUE)+scale_colour_manual(values=c('TRUE'='#7560A8','FALSE'=grey),guide='none')+scale_size_continuous(range=c(1.4,3.6),name='Prior-export log2FC')+labs(title='Prior selection of downstream markers',subtitle='Existing ranking retained; not therapeutic efficacy or independent target validation',x='External depletion (-median log2FC)',y='Platform depletion (-log2FC)')
savep(p,'Supplementary_Figure_6_downstream_ranking',135)
d<-old('SF7_effects');gs<-sort(unique(d$gene));d$context<-factor(d$context,levels=c('Organoid TA',families));d$gene<-factor(d$gene,levels=rev(gs));d$mark<-ifelse(is.na(d$effect),'NE',ifelse(d$q<.05,'*',''));half<-ceiling(length(gs)/2)
parts<-lapply(split(gs,ceiling(seq_along(gs)/half)),function(g){ggplot(d[d$gene%in%g,],aes(context,gene,fill=effect))+geom_tile(colour='white')+geom_text(aes(label=mark),size=2.4)+nes_scale(d$effect,'log2FC')+scale_x_discrete(labels=c('Prior export TA',flabs))+labs(x=NULL,y=NULL)+theme(axis.text.x=element_text(angle=50,hjust=1),legend.position='bottom')})
savep(compose(ggdraw()+draw_label('Complete selected core: prior export and external gene effects',size=9,fontface='bold'),plot_grid(plotlist=parts,ncol=2,labels=c('a','b'),label_size=10),ggdraw()+draw_label('Alphabetical genes; * each source gene-level FDR < 0.05. External data informed selection.',size=7),rh=c(.05,.90,.05)),'Supplementary_Figure_7_complete_core_effects',230)
tn<-old('F4B');tn$module<-factor(tn$module,levels=rev(modids));p1<-ggplot(tn,aes(NES,module,shape=screen_padj<.05))+geom_vline(xintercept=0,colour=grey,linetype=2)+geom_point(size=2,colour=blue)+scale_y_discrete(labels=modlab)+scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),guide='none')+labs(title='Prior export-based TNFSF12',subtitle='Filled: original screen q < 0.05',x='Raw NES',y=NULL)
ef<-old('F4C');p2<-ggplot(ef,aes(effect_size,module,colour=celltype,shape=padj<.05))+scale_colour_manual(values=pal_family,labels=ctlab)+geom_vline(xintercept=0,colour=grey,linetype=2)+geom_point(position=position_dodge(width=.45,orientation='y'),size=2)+scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),guide='none')+labs(title='Published programme context',subtitle='Filled: publisher q < 0.05',x="Published Cohen's d",y=NULL,colour='Family')+theme(legend.position='bottom')+guides(colour=guide_legend(ncol=2))
ref<-old('SF14C_REFERENCE');ref$ct<-factor(ref$ct,levels=names(pal_family));ref$module<-factor(ref$module,levels=rev(modids));p3<-ggplot(ref,aes(ct,module,fill=NES))+geom_tile(colour='white')+geom_text(aes(label=mark),size=2.6)+nes_scale(ref$NES)+scale_x_discrete(labels=ctlab)+scale_y_discrete(labels=modlab)+labs(title='Count-based assay reference without cytomix',subtitle='Raw NES; * separate reference 16-test FDR < 0.05',x=NULL,y=NULL)
savep(compose(compose(p1,p2,labs=c('a','b'),ncol=2),p3,labs=c('','c'),rh=c(.58,.42)),'Supplementary_Figure_14_original_and_assay_context',180)

cat('SUPPLEMENT_READER_PLOTS_PASS\n')
