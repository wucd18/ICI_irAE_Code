#!/usr/bin/env Rscript
# Only reads the exact TARGETED_FINAL snapshot. One authorized display lm in F6B.

# Recorded display annotations are external scientific metadata, never code constants.
archived_label <- function(key) {
 p<-Sys.getenv('ICI_HISTORICAL_DISPLAY_LABELS',unset='')
 if(!nzchar(p)||!file.exists(p))stop('BLOCKED_INPUT: ICI_HISTORICAL_DISPLAY_LABELS')
 labels<-jsonlite::fromJSON(p)
 if(is.null(labels[[key]])||length(labels[[key]])!=1)stop(paste('Missing archived label',key))
 as.character(labels[[key]])
}

options(warn=1)
args<-commandArgs(TRUE);stopifnot(length(args)==1)
O<-normalizePath(args[1],winslash='/');I<-file.path(O,'visual_input_snapshot')
source(file.path(O,'presentation_code/theme_final.R'))
source(file.path(O,'presentation_code/csv_io.R'))
required_archive_path<-function(key){p<-Sys.getenv(key,unset='');if(!nzchar(p)||!dir.exists(p))stop(paste('BLOCKED_INPUT',key));normalizePath(p,winslash='/')}
rd<-function(p)read.csv(file.path(I,p),check.names=FALSE,na.strings=c('NA',''),stringsAsFactors=FALSE)
fp<-function(p)read.csv(file.path(O,'source_data/final_panels',paste0(p,'.csv')),check.names=FALSE,na.strings=c('NA',''),stringsAsFactors=FALSE)
old<-function(p)rd(paste0('source_data/panel_sources/',p,'.csv'))
tab<-function(n) {p<-list.files(file.path(I,'supplement/full_source_tables'),pattern=sprintf('^Source_%03d_',n));stopifnot(length(p)==1);rd(paste0('supplement/full_source_tables/',p))}
writecsv<-function(d,p){d<-as.data.frame(d);for(k in names(d))if(is.list(d[[k]]))d[[k]]<-vapply(d[[k]],function(x)as.character(jsonlite::toJSON(x,auto_unbox=TRUE)),character(1));write.csv(d,file.path(O,p),row.names=FALSE,na='NA')}
for(n in c('figures','review/plot_build','source_data/display','logs'))dir.create(file.path(O,n),recursive=TRUE,showWarnings=FALSE)
theme_update(legend.text=element_text(size=7),strip.text=element_text(size=7.5),plot.margin=margin(5,6,5,6))
flabs<-c('Absorptive','Immature','Mature absorptive','Secretory','GSE210037\nsample pseudobulk')
modlabs<-c('IFN reference','Disease-up','IFN-gene-excluded\nresidual','Disease-depleted')
modlab<-function(x)unname(setNames(modlabs,modids)[x])
dims<-list();checks<-list()
savep<-function(p,id,h){
 base<-file.path(O,'figures',id)
 ggsave(paste0(base,'.pdf'),p,width=180,height=h,units='mm',device=cairo_pdf,bg='white')
 ggsave(paste0(base,'.tiff'),p,width=180,height=h,units='mm',dpi=600,device=grDevices::tiff,compression='lzw',bg='white')
 ggsave(paste0(base,'_preview.png'),p,width=180,height=h,units='mm',dpi=220,bg='white')
 dims[[id]]<<-data.frame(figure=id,width_mm=180,height_mm=h,tiff_dpi=600)
 message('Rendered ',id)
}
audit<-function(id,p,d,layer,value,xcol='x',transform=identity){
 b<-ggplot_build(p);saveRDS(b,file.path(O,'review/plot_build',paste0(id,'.rds')))
 for(k in seq_along(b$data))writecsv(b$data[[k]],paste0('review/plot_build/',id,'_layer_',k,'.csv'))
 got<-sort(b$data[[layer]][[xcol]]);want<-sort(transform(d[[value]]))
 stopifnot(length(got)==length(want),isTRUE(all.equal(got,want,tolerance=1e-12,check.attributes=FALSE)))
 checks[[id]]<<-data.frame(panel=id,objects=nrow(d),field=value,max_abs_difference=max(abs(got-want)),status='PASS')
 writecsv(d,paste0('source_data/display/',id,'.csv'))
}
whitegap_y<-function(x) {z<-match(x,rev(families));ifelse(z==1,.55,z)}
# Figure 1: editable R paths and real coverage, no fictional observations.
h<-fp('F2A');v<-fp('F3A')
cov<-bind_rows(h%>%distinct(dataset,organ,celltype,n_case,n_control)%>%mutate(group=ifelse(organ=='heart','Heart','Colon')),
 v%>%distinct(validation_dataset,celltype,n_case,n_control)%>%rename(dataset=validation_dataset)%>%mutate(group=ifelse(dataset=='GSE210037','GSE210037','GSE206300')))
cov$label<-c('Endothelial','Fibroblast','Mural','Endothelial','Epithelial','Stromal',flabs[1:4],'Sample pseudobulk')
# Explicit ordered join uses source identities, never hand-entered observed counts.
heart_order<-c('Endothelial cells','Fibroblasts','Mural cells');colon_order<-c('Endothelial','Epithelial','Mesenchymal Stromal')
cov<-bind_rows(cov%>%filter(group=='Heart')%>%arrange(match(celltype,heart_order)),cov%>%filter(group=='Colon')%>%arrange(match(celltype,colon_order)),cov%>%filter(group=='GSE206300')%>%arrange(match(celltype,families)),cov%>%filter(group=='GSE210037'))
cov$label<-c('Endothelial','Fibroblast','Mural','Endothelial','Epithelial','Stromal',flabs[1:4],'Sample pseudobulk')
cov$y<-rev(seq_len(nrow(cov)))+c(rep(1.2,3),rep(.8,3),rep(.4,4),0)
cv<-cov%>%pivot_longer(c(n_case,n_control),names_to='condition',values_to='n')%>%mutate(yy=y+ifelse(condition=='n_case',.14,-.14))
writecsv(cov,'source_data/display/F1B_coverage.csv')
ng<-read.csv(file.path(required_archive_path('ICI_RECORDED_ORGANOID_RESULTS'),'module_GSEA.csv'),check.names=FALSE)
inv<-read.csv(file.path(required_archive_path('ICI_RECORDED_ORGANOID_RESULTS'),'model_inventory.csv'),check.names=FALSE)
md<-rd('organoid_donor_extension/source_data/pseudobulk_metadata.csv')
donors<-sort(unique(md$donor));model_plan<-nrow(unique(ng[c('ligand','celltype')]))
counts<-data.frame(kind=c('Count models','Set comparisons'),estimable=c(nrow(inv),sum(ng$status=='ESTIMATED')),planned=c(model_plan,nrow(ng)))
writecsv(counts,'source_data/display/F1C_coverage.csv')
writecsv(data.frame(donor=donors),'source_data/display/F1C_donors.csv')
# Arrows are workflow arrows, placed behind the vector icons.
p1a<-ggplot()+annotate('segment',x=c(19,44,69),xend=c(27,52,77),y=63,yend=63,linewidth=.35,colour=grey,arrow=arrow(length=unit(1.5,'mm')))+
 annotate('segment',x=12,xend=c(12,36),y=37,yend=c(15,15),linetype=2,colour=grey,linewidth=.3)+
 annotate('text',x=c(12,36),y=5,label=c('TCR support','Lung immune sensitivity'),size=2.5,colour=ink)
theta<-seq(0,2*pi,length.out=150)
heart<-data.frame(x=9+.20*16*sin(theta)^3,y=66+.65*(13*cos(theta)-5*cos(2*theta)-2*cos(3*theta)-cos(4*theta)))
colon<-data.frame(x=c(15,14,14,15,18,19,19,18,16,16),y=c(57,59,70,72,72,70,60,58,58,55))
p1a<-p1a+geom_polygon(data=heart,aes(x,y),fill='#EDD0CE',colour=red,linewidth=.35)+geom_path(data=colon,aes(x,y),colour=blue,linewidth=1.1,lineend='round',linejoin='round')+
 annotate('rect',xmin=30,xmax=42,ymin=53,ymax=74,fill='white',colour=grey,linewidth=.3)+
 annotate('segment',x=32,xend=c(39,40,37,39),y=c(70,66,62,58),yend=c(70,66,62,58),colour=unname(pal_module),linewidth=1.4)+
 annotate('rect',xmin=55,xmax=67,ymin=54,ymax=73,fill='white',colour=grey,linewidth=.4)+
 geom_point(data=expand.grid(x=seq(57,65,length.out=3),y=seq(57,70,length.out=4)),aes(x,y),shape=21,fill='#D8E8E6',colour='#2E8079',size=1.6,stroke=.3)+
 annotate('segment',x=80,xend=80,y=54,yend=75,colour=ink,linewidth=.3)+annotate('segment',x=80,xend=93,y=54,yend=54,colour=ink,linewidth=.3)+
 geom_path(data=data.frame(x=c(81,84,87,90,93),y=c(72,72,66,66,60)),aes(x,y),colour='#7560A8',linewidth=.6)+
 annotate('text',x=c(13,36,61,86),y=91,label=c('Human tissues','Frozen programmes','Paired organoids','Clinical associations'),size=3,fontface='bold')+
 annotate('text',x=c(13,36,61,86),y=42,label=c('Heart and colon','Independent transfer','Donor-level counts','Tumour guardrails'),size=2.5)+
 coord_cartesian(xlim=c(0,100),ylim=c(0,100),clip='off')+theme_void()+theme(plot.margin=margin(6,5,5,10))
p1b<-ggplot(cv,aes(n,yy,colour=condition))+geom_segment(aes(x=0,xend=n,yend=yy),linewidth=.5,alpha=.55)+geom_point(size=1.7)+geom_text(aes(label=n),nudge_x=.3,hjust=0,size=2.5,show.legend=FALSE)+
 scale_colour_manual(values=c(n_case=red,n_control=grey),labels=c('ICI cases','Controls'),name=NULL)+scale_y_continuous(breaks=cov$y,labels=paste0(cov$group,': ',cov$label),expand=expansion(add=.5))+
 scale_x_continuous(limits=c(0,max(cv$n)+2),breaks=pretty(c(0,max(cv$n)),4))+labs(title='Biological replicates in core comparisons',x='Patients or released tissue samples',y=NULL)+theme(legend.position='bottom',axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=7))
donorplot<-ggplot(data.frame(x=seq_along(donors),donor=donors),aes(x,1))+
 geom_point(shape=21,size=3,fill='#D8E8E6',colour='#2E8079')+geom_segment(aes(xend=x,y=.78,yend=.45),linewidth=2,colour='#2E8079')+
 geom_text(aes(label=donor),y=.16,size=2.5)+coord_cartesian(xlim=c(.5,length(donors)+.5),ylim=c(0,1.3))+theme_void()+labs(title=paste(length(donors),'paired human donors'))+theme(plot.title=element_text(size=8.5,face='bold'),plot.margin=margin(5,5,5,20))
bars<-ggplot(counts,aes(estimable,factor(kind,levels=rev(kind))))+geom_col(aes(x=planned),fill='#E4E8EB',width=.46)+geom_col(fill=blue,width=.46)+
 geom_text(aes(label=paste0(estimable,' / ',planned)),hjust=-.12,size=2.5)+scale_x_continuous(limits=c(0,max(counts$planned)*1.35),breaks=pretty(c(0,max(counts$planned)),3))+
 labs(x='Estimable / planned comparisons',y=NULL,subtitle='Blue: estimable; pale grey: unavailable')+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.subtitle=element_text(size=7))
p1c<-compose(donorplot,bars,ggdraw()+draw_label('Coverage counts are not independent donors.\nEstimability does not imply a positive result.',size=7,x=0,hjust=0),rh=c(.32,.5,.18))
savep(compose(p1a,compose(p1b,p1c,labs=c('b','c'),ncol=2,rw=c(.63,.37)),labs=c('a',''),rh=c(.37,.63)),'Figure_1_study_design',145)
# Figure 2: unchanged numerical coordinates and controlled white gaps.
s<-fp('F2B');l<-fp('F2C');stats<-old('F2B_statistics')
shared<-s$pathway[s$class=='Shared up'];divergent<-s$pathway[s$class=='Context-divergent'];path_order<-c(shared,divergent)
h$x<-match(h$compartment,unique(h$compartment));h$x<-h$x+ifelse(h$x>3,.28,0)
h$y<-match(h$pathway,rev(path_order));h$y<-h$y+ifelse(h$y>length(divergent),.5,0)
hy<-sort(unique(h$y));h$mark[is.na(h$mark)&!is.na(h$NES)]<-''
p2a<-ggplot(h,aes(x,y,fill=NES))+geom_tile(width=.98,height=.96,colour=NA)+geom_text(aes(label=mark),size=2.6)+nes_scale(h$NES)+
 scale_y_continuous(breaks=hy,labels=pathlab(rev(path_order)),expand=expansion(add=.7))+
 scale_x_continuous(breaks=sort(unique(h$x)),labels=c('Endothelial','Fibroblast','Mural','Endothelial','Epithelial','Stromal'),position='bottom',expand=expansion(add=.5))+
 annotate('text',x=c(2,5.28),y=max(h$y)+1.25,label=c('Heart','Colon'),size=2.8,fontface='bold')+
 labs(title='Shared inflammatory visibility and receiver-cell remodelling',x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank())
s$lab<-ifelse(s$class=='Shared up',gsub(' response| / NF-kB| / JAK / STAT3| / STAT5','',pathlab(s$pathway)),ifelse(s$class=='Context-divergent',as.character(match(s$pathway,divergent)),''))
writecsv(data.frame(index=seq_along(divergent),pathway=divergent,label=pathlab(divergent)),'source_data/display/F2B_label_key.csv')
xy<-range(c(s$heart_median_NES,s$colon_median_NES));lims<-xy+c(-.10,.40)*diff(xy)
edge<-s[s$class=='Shared up',];edge<-edge[order(edge$colon_median_NES,decreasing=TRUE),];edge$tx<-xy[2]+diff(xy)*.14;edge$ty<-seq(xy[2]+.15,.6,length.out=nrow(edge));edge$short<-c('IFN-gamma','IFN-alpha','G2M','Allograft','Inflamm.','Complement')[match(edge$pathway,c('HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_G2M_CHECKPOINT','HALLMARK_ALLOGRAFT_REJECTION','HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_COMPLEMENT'))]
s$lab[s$class=='Shared up']<-''
p2b<-ggplot(s,aes(heart_median_NES,colon_median_NES,colour=class))+geom_hline(yintercept=0,colour=grey,linewidth=.2)+geom_vline(xintercept=0,colour=grey,linewidth=.2)+geom_point(size=1.65)+
 geom_text_repel(aes(label=lab),size=2.5,seed=1,max.overlaps=Inf,box.padding=.25,point.padding=.2,force=2,max.time=3,max.iter=20000,min.segment.length=0,segment.size=.18,show.legend=FALSE)+
 geom_segment(data=edge,aes(xend=tx-.08,yend=ty),linewidth=.18,show.legend=FALSE)+geom_text(data=edge,aes(x=tx,y=ty,label=short),hjust=0,size=2.5,show.legend=FALSE)+
 scale_colour_manual(values=c('Shared up'=red,'Context-divergent'=blue,Other=grey))+coord_equal(xlim=lims,ylim=lims,clip='off')+labs(title='All Hallmark programmes',subtitle=paste0('rho = ',sprintf('%.4f',stats$spearman_rho),'; P = ',fmt(stats$p)),x='Heart median NES',y='Colon median NES',colour=NULL)+theme(legend.position='bottom')+guides(colour=guide_legend(ncol=2))
l$module<-factor(l$pathway,levels=rev(modids));l$cell<-factor(l$celltype);l$y<-as.numeric(l$module)+seq(-.25,.25,length.out=nlevels(l$cell))[l$cell]
p2c<-ggplot(l,aes(NES,y,colour=celltype,shape=padj<.05))+geom_vline(xintercept=0,linetype=2,colour=grey,linewidth=.25)+geom_point(size=1.9)+scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),guide='none')+
 scale_y_continuous(breaks=1:4,labels=rev(c('IFN (+)','Disease-up (+)','Residual (+)','Depleted (-)')))+scale_colour_manual(values=c('#2B6E98','#B84B46','#7560A8','#B58B20'),labels=c('Alveolar macro.','CD4 T','Mono-derived macro.','cDCs'))+
 labs(title='Lung immune sensitivity',subtitle='Source-confounded; filled: FDR < 0.05',x='NES (CIP vs external healthy BALF)',y=NULL,colour=NULL)+theme(legend.position='bottom',axis.title.x=element_text(size=7))+guides(colour=guide_legend(ncol=2))
savep(compose(p2a,compose(p2b,p2c,labs=c('b','c'),ncol=2,rw=c(1.03,.97)),labs=c('a',''),rh=c(.53,.47)),'Figure_2_receiver_programmes',175)
audit('F2B',p2b,s,3,'heart_median_NES');audit('F2C',p2c,l,2,'NES');writecsv(h,'source_data/display/F2A.csv')
# Figure 3: dot matrix, original AUC CIs, actual membership partition.
a<-fp('F3B');v$y<-whitegap_y(v$validation);a$y<-whitegap_y(a$validation);v$module<-factor(v$module,levels=modids);a$module<-factor(a$module,levels=modids)
fy<-whitegap_y(families)
p3a<-ggplot(v,aes(module,y,fill=oriented_NES))+geom_point(shape=21,size=4.1,stroke=.3,colour=ink)+geom_text(aes(label=mark),size=2.5)+nes_scale(v$oriented_NES,'Oriented NES')+
 scale_x_discrete(labels=modlab)+scale_y_continuous(breaks=fy,labels=flabs,expand=expansion(add=.65))+labs(title='Frozen programmes across independent measurements',x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank())
p3b<-ggplot(a,aes(auc_disease_higher,y,colour=module))+geom_vline(xintercept=.5,linetype=2,colour=grey,linewidth=.25)+geom_segment(aes(x=bootstrap_auc_ci_low,xend=bootstrap_auc_ci_high,yend=y),linewidth=.4)+geom_point(size=1.8)+
 facet_wrap(~module,ncol=2,labeller=as_labeller(setNames(modlabs,modids)))+scale_colour_manual(values=unname(pal_module),guide='none')+scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1),expand=expansion(mult=c(.035,.035)))+
 scale_y_continuous(breaks=fy,labels=flabs,expand=expansion(add=.65))+labs(x='AUC',y=NULL)+theme(panel.spacing.x=unit(7,'mm'),panel.spacing.y=unit(4,'mm'))
m<-fp('F3C');sets<-lapply(split(m$gene,m$module),unique);overlap<-intersect(sets[[modids[2]]],sets[[modids[1]]]);stopifnot(setequal(setdiff(sets[[modids[2]]],overlap),sets[[modids[3]]]))
sz<-vapply(sets,length,integer(1));part<-data.frame(segment=c('Residual','IFN overlap'),start=c(0,sz[[modids[3]]]),end=c(sz[[modids[3]]],sz[[modids[2]]]))
writecsv(data.frame(module=names(sz),members=sz),'source_data/display/F3C_sizes.csv');writecsv(data.frame(gene=overlap),'source_data/display/F3C_excluded.csv')
p3c<-ggplot(part)+geom_rect(aes(xmin=start,xmax=end,ymin=.8,ymax=1.2,fill=segment),colour='white',linewidth=.2)+scale_fill_manual(values=c(Residual='#7560A8','IFN overlap'='#A44D79'),guide='none')+
 annotate('text',x=0,y=1.55,hjust=0,label=paste0('Disease-up: ',sz[[modids[2]]],' members'),size=2.7,fontface='bold')+
 annotate('text',x=sz[[modids[3]]]/2,y=1,label=paste0(sz[[modids[3]]],' residual (',sprintf('%.2f',100*sz[[modids[3]]]/sz[[modids[2]]]),'%)'),size=2.6,colour='white')+
 annotate('text',x=sz[[modids[2]]]+3,y=1,hjust=0,label=paste0(length(overlap),' IFN overlap'),size=2.5)+
 annotate('text',x=0,y=.3,hjust=0,label=paste0('Separate reference sizes: IFN ',sz[[modids[1]]],'   |   Disease-depleted ',sz[[modids[4]]]),size=2.5)+
 annotate('text',x=0,y=-.12,hjust=0,label='Excluding IFN-set genes does not establish IFN independence',size=2.5)+
 coord_cartesian(xlim=c(0,sz[[modids[2]]]*1.5),ylim=c(-.3,1.8),clip='off')+labs(title='Definition of the IFN-gene-excluded set')+theme_void()+theme(plot.title=element_text(size=8.5,face='bold'),plot.margin=margin(8,12,4,30))
savep(compose(p3a,p3b,p3c,labs=c('a','b','c'),rh=c(.29,.49,.22)),'Figure_3_frozen_module_transfer',175)
audit('F3A',p3a,v,1,'y','y');audit('F3B',p3b,a,3,'auc_disease_higher')
# Figure 4: fixed set lollipops and donor changes remain separate metrics.
n<-fp('F4A');n$module<-factor(n$module,levels=modids);n$celltype<-factor(n$celltype,levels=rev(names(pal_family)))
p4a<-ggplot(n,aes(NES,celltype,colour=celltype))+geom_vline(xintercept=0,linetype=2,colour=grey,linewidth=.25)+geom_segment(aes(x=0,xend=NES,yend=celltype),linewidth=.35)+geom_point(aes(shape=FDR_predeclared_family<.05),size=2)+
 facet_wrap(~module,ncol=2,labeller=as_labeller(setNames(modlabs,modids)))+scale_colour_manual(values=pal_family,guide='none')+scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),labels=c('FDR >= 0.05','FDR < 0.05'),name=NULL)+scale_y_discrete(labels=ctlab)+
 labs(title='TNFSF12 responses in fixed ICI programmes',subtitle=paste0(unique(n$n_donors),' paired donors; ',nrow(n),'-test family'),x='NES',y=NULL)+theme(legend.position='bottom',panel.spacing.x=unit(3,'mm'),panel.spacing.y=unit(4,'mm'),plot.title=element_text(size=8))
d<-fp('F4B');d$celltype<-factor(d$celltype,levels=rev(names(pal_family)));d$donor<-factor(d$donor,levels=names(pal_donor));d$y<-as.numeric(d$celltype)+seq(-.18,.18,length.out=nlevels(d$donor))[d$donor]
means<-d%>%group_by(celltype)%>%summarise(mean=mean(change),.groups='drop')
p4b<-ggplot(d,aes(change,y,colour=donor))+geom_vline(xintercept=0,linetype=2,colour=grey,linewidth=.25)+geom_segment(aes(x=0,xend=change,yend=y),linewidth=.3)+geom_point(size=2)+
 geom_segment(data=means,aes(x=mean,xend=mean,y=as.numeric(celltype)-.1,yend=as.numeric(celltype)+.1),inherit.aes=FALSE,linewidth=.5,colour=ink)+scale_colour_manual(values=pal_donor)+scale_y_continuous(breaks=1:4,labels=rev(ctlab),expand=expansion(add=.4))+
 labs(title='Disease-depleted programme',subtitle='Donor changes; black tick = mean',x='Paired change in\nmean expression rank',y=NULL,colour='Donor')+theme(legend.position='bottom',plot.title=element_text(size=8))+guides(colour=guide_legend(ncol=1))
rc<-fp('F4C');rc$x<-match(rc$context,families);rc$x<-rc$x+ifelse(rc$x==5,.35,0);rc$signature<-factor(rc$pathway,levels=rev(c('TNFSF12_ORGANOID_UP','TNFSF12_ORGANOID_DOWN','TNFSF12_FUNCTION_RESCUE_CORE')))
p4c<-ggplot(rc,aes(x,signature,fill=NES))+geom_tile(width=.98,height=.96,colour=NA)+geom_text(aes(label=mark,colour=abs(NES)>max(abs(NES))*.63),size=2.7)+scale_colour_manual(values=c('FALSE'=ink,'TRUE'='white'),guide='none')+nes_scale(rc$NES)+
 scale_x_continuous(breaks=c(1:4,5.35),labels=flabs,expand=expansion(add=.5))+scale_y_discrete(labels=c(archived_label('selected_core_label'),archived_label('selected_down_label'),archived_label('selected_up_label')))+labs(title='Previously selected sets in external disease tissue',subtitle='GSE206300 families and GSE210037 whole-tissue samples; * source FDR < 0.05',x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank())
savep(compose(compose(p4a,p4b,labs=c('a','b'),ncol=2,rw=c(.62,.38)),p4c,labs=c('','c'),rh=c(.66,.34)),'Figure_4_organoid_reciprocity',180)
audit('F4A',p4a,n,3,'NES');audit('F4B',p4b,d,3,'change');writecsv(rc,'source_data/display/F4C.csv')
# Figure 5: aligned clinical estimates and independent FDR column for SORL1.
r<-fp('F5A');sv<-fp('F5B');ord<-c(r$feature_id[r$category=='GENE'],r$feature_id[r$category=='AXIS']);pos<-setNames(rev(seq_along(ord)),ord);pos<-pos+ifelse(startsWith(names(pos),'GENE'),.8,0)
r$y<-unname(pos[r$feature_id]);sv$y<-unname(pos[sv$feature_id]);featurelab<-sub('^GENE_|^AXIS_','',ord)
ylab<-paste0(ifelse(startsWith(ord,'GENE'),'Gene  ','Axis  '),featurelab)
forest_theme<-theme(axis.ticks.y=element_blank(),axis.line.y=element_blank(),plot.margin=margin(4,10,4,4))
p5a<-ggplot(r,aes(estimate,y))+geom_vline(xintercept=0,linetype=2,colour=grey,linewidth=.25)+geom_segment(aes(x=ci_low,xend=ci_high,yend=y),colour=blue,linewidth=.35)+geom_point(aes(shape=fdr<.05),size=1.8,colour=blue)+
 scale_y_continuous(breaks=unname(pos),labels=ylab,expand=expansion(add=.6))+scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),guide='none')+labs(title='Tumour response',subtitle='Filled: original FDR < 0.05',x="Hedges' g (95% CI)",y=NULL)+forest_theme
p5b<-ggplot(sv,aes(hazard_ratio,y))+geom_vline(xintercept=1,linetype=2,colour=grey,linewidth=.25)+geom_segment(aes(x=ci_low_hr,xend=ci_high_hr,yend=y),colour='#7560A8',linewidth=.35)+geom_point(aes(shape=fdr<.05),size=1.8,colour='#7560A8')+
 scale_y_continuous(breaks=unname(pos),labels=ylab,expand=expansion(add=.6))+scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),guide='none')+scale_x_log10(breaks=c(.6,.8,1,1.25,1.5))+labs(title='Overall survival',subtitle='Same feature order',x='Hazard ratio (95% CI, log scale)',y=NULL)+forest_theme
so<-fp('F5C');so$y<-rev(seq_len(nrow(so)))+ifelse(so$source_group=='Donor-count organoids',.8,ifelse(so$source_group=='External patient families',.4,0))
src_cols<-c('Donor-count organoids'=blue,'External patient families'='#7560A8','Separate platform'='#2E8079')
ctx<-so$context;ctx<-gsub('Transit_Amplifying','TA',ctx);ctx<-gsub(' epithelial','',ctx);ctx[so$source_group=='Separate platform']<-'Sample pseudobulk'
yr<-range(so$y)+c(-.5,.5)
grouptext<-data.frame(y=c(mean(so$y[so$source_group=='Donor-count organoids']),mean(so$y[so$source_group=='External patient families']),mean(so$y[so$source_group=='Separate platform'])),label=c('Donor-count\norganoids','GSE206300\nfour families','GSE210037'))
left<-ggplot(grouptext,aes(0,y))+geom_text(aes(label=label),hjust=0,size=2.5,fontface='bold')+scale_y_continuous(limits=yr,expand=c(0,0))+coord_cartesian(xlim=c(0,1),clip='off')+theme_void()+theme(plot.margin=margin(20,0,25,0))
effect<-ggplot(so,aes(estimate,y,colour=source_group))+geom_vline(xintercept=0,linetype=2,colour=grey,linewidth=.25)+geom_segment(aes(x=0,xend=estimate,yend=y),linewidth=.4)+geom_point(aes(shape=q<.05),size=2)+
 scale_colour_manual(values=src_cols,guide='none')+scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),guide='none')+scale_y_continuous(breaks=so$y,labels=ctx,limits=yr,expand=c(0,0))+labs(title='SORL1 measurements',x='log2 fold change',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.margin=margin(5,5,5,3))
qcol<-ggplot(so,aes(0,y))+geom_text(aes(label=fmt(q)),hjust=0,size=2.5)+scale_y_continuous(limits=yr,expand=c(0,0))+coord_cartesian(xlim=c(0,1),clip='off')+labs(title='Gene FDR')+theme_void()+theme(plot.title=element_text(size=8,face='bold'),plot.margin=margin(5,5,25,0))
aligned<-align_plots(left,effect,qcol,align='h',axis='tb');p5c<-plot_grid(plotlist=aligned,nrow=1,rel_widths=c(.18,.7,.12))
savep(compose(compose(p5a,p5b,labs=c('a','b'),ncol=2),p5c,labs=c('','c'),rh=c(.58,.42)),'Figure_5_tumour_guardrails',175)
audit('F5A',p5a,r,3,'estimate');audit('F5B',p5b,sv,3,'hazard_ratio',transform=log10);audit('F5C',effect,so,3,'estimate')
# Figure 6: original patient data and tests; explicit display-only OLS.
t<-fp('F6A');t$group<-factor(t$group,levels=c('Healthy','ICI colitis'));t$metric<-factor(t$metric,levels=c('expanded_cell_fraction_ge2','normalized_clonality','cross_cd4_cd8_cell_fraction'))
ts<-tab(35)%>%filter(analysis=='colon_CPI_colitis_vs_HC_patient_level',metric%in%levels(t$metric));ts$metric<-factor(ts$metric,levels=levels(t$metric))
tlabels<-c(expanded_cell_fraction_ge2='Expanded-cell fraction',normalized_clonality='Normalized clonality',cross_cd4_cd8_cell_fraction='CD4-CD8 shared fraction')
p6a<-ggplot(t,aes(group,value,fill=group))+geom_boxplot(width=.5,outlier.shape=NA,linewidth=.3,alpha=.4)+geom_point(position=position_jitter(width=.08,height=0,seed=1),shape=21,colour=ink,size=1.8,stroke=.25)+
 geom_text(data=ts,aes(x=1.5,y=Inf,label=paste0('P=',fmt(p),'  FDR=',fmt(fdr_within_analysis))),inherit.aes=FALSE,vjust=1,size=2.5)+facet_wrap(~metric,nrow=1,scales='free_y',labeller=as_labeller(tlabels))+
 scale_fill_manual(values=pal_group,guide='none')+scale_y_continuous(expand=expansion(mult=c(.05,.23)))+labs(title='Patient-level TCR expansion provides supportive immune-pressure evidence',x=NULL,y='Observed patient metric')+theme(axis.text.x=element_text(angle=18,hjust=1))
tc<-fp('F6B');cs<-old('F6B_stats');stopifnot(!anyDuplicated(tc$patient_id),all(tc$disease=='CPI_colitis'),setequal(tc$patient_id,unique(t$patient_id[t$group=='ICI colitis'])))
stopifnot(all(is.finite(tc$cross_cd4_cd8_cell_fraction)),all(is.finite(tc$ICI_COLITIS_EPITHELIAL_LOSS_STRICT)))
# The only model fit in this entry: no observations excluded; no new test reported.
display_fit<-lm(ICI_COLITIS_EPITHELIAL_LOSS_STRICT~cross_cd4_cd8_cell_fraction,data=tc)
prediction<-data.frame(cross_cd4_cd8_cell_fraction=seq(min(tc$cross_cd4_cd8_cell_fraction),max(tc$cross_cd4_cd8_cell_fraction),length.out=160))
prediction<-cbind(prediction,as.data.frame(predict(display_fit,newdata=prediction,interval='confidence',level=.95)))
writecsv(prediction,'source_data/display/F6B_OLS_mean_band.csv')
jsonlite::write_json(list(purpose='Descriptive display only',patient_ids=tc$patient_id,n=nrow(tc),formula='unreversed mean expression rank ~ CD4-CD8 shared-clone cell fraction',band='pointwise 95% confidence interval for conditional OLS mean response; not a Spearman interval or prediction interval',coefficients=coef(display_fit),residual_df=df.residual(display_fit),original_spearman=cs),file.path(O,'review/F6B_display_fit.json'),pretty=TRUE,auto_unbox=TRUE)
p6b<-ggplot(tc,aes(cross_cd4_cd8_cell_fraction,ICI_COLITIS_EPITHELIAL_LOSS_STRICT))+
 geom_ribbon(data=prediction,aes(x=cross_cd4_cd8_cell_fraction,ymin=lwr,ymax=upr),inherit.aes=FALSE,fill='#DDE4E8')+geom_line(data=prediction,aes(x=cross_cd4_cd8_cell_fraction,y=fit),inherit.aes=FALSE,colour='#9AA4AD',linewidth=.55)+
 geom_point(colour='#356B91',size=2.1)+geom_text_repel(aes(label=patient_id),seed=1,size=2.5,max.overlaps=Inf,min.segment.length=0,segment.size=.2)+
 labs(title='Exploratory epithelial association',subtitle=paste0('Spearman rho=',fmt(cs$spearman_rho),'; P=',fmt(cs$p),'\nOriginal 16-test FDR=',fmt(cs$fdr)),x='CD4-CD8 shared-clone cell fraction',y='Unreversed mean expression rank')
rare<-fp('F6C');rare$label<-factor(rare$label,levels=rev(rare$label))
p6c<-ggplot(rare,aes(median_difference,label))+geom_vline(xintercept=0,colour=grey,linetype=2,linewidth=.25)+geom_segment(aes(x=0,xend=median_difference,yend=label),colour='#B58B20',linewidth=.55)+geom_point(size=2.7,colour='#B58B20')+
 geom_text(aes(x=0,label=paste0('P=',fmt(p),'  FDR=',fmt(fdr_within_analysis))),vjust=-1.5,hjust=0,size=2.5)+scale_y_discrete(labels=function(x)gsub('Rarefied ','Rarefied\n',x),expand=expansion(add=.9))+scale_x_continuous(expand=expansion(mult=c(.06,.2)))+
 geom_text(aes(x=0,label=paste0('Median difference = ',sprintf('%.5f',median_difference))),vjust=2,hjust=0,size=2.5)+geom_text(aes(x=0,label=paste0('Median difference = ',sprintf('%.5f',median_difference))),vjust=2,hjust=0,size=2.5)+labs(title='Depth-controlled sensitivity',subtitle='Each effect retains its own metric units',x='Median difference\n(ICI colitis - healthy)',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
savep(compose(p6a,compose(p6b,p6c,labs=c('b','c'),ncol=2,rw=c(.59,.41)),labs=c('a',''),rh=c(.43,.57)),'Figure_6_TCR_support',165)
audit('F6A',p6a,t,2,'value','y');audit('F6B',p6b,tc,3,'cross_cd4_cd8_cell_fraction');audit('F6C',p6c,rare,3,'median_difference')
writecsv(bind_rows(dims),'review/main_figure_dimensions.csv');writecsv(bind_rows(checks),'review/plot_coordinate_checks.csv')
for(id in c('F1A','F1B','F1C','F2A','F3C','F4C')) {p<-get(tolower(sub('F','p',id)));saveRDS(p,file.path(O,'review/plot_build',paste0(id,'_plot.rds')))}
capture.output(sessionInfo(),file=file.path(O,'logs/presentation_R_sessionInfo.txt'))
cat('VISUAL_MAIN_PLOTS_PASS; only authorized F6B display OLS was fitted\n')
