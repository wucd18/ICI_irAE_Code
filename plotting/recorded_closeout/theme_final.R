suppressPackageStartupMessages({library(ggplot2);library(cowplot);library(dplyr);library(tidyr);library(ggrepel)})
if(.Platform$OS.type=='windows')grDevices::windowsFonts(Arial=grDevices::windowsFont('Arial'))
pal_group<-c(Healthy='#9AA4AD','ICI colitis'='#B84B46',UC='#2E8079')
pal_module<-c(IFN='#A44D79',Up='#B58B20',Residual='#7560A8',Depleted='#2B6E98')
pal_family<-c(Stem='#2B6E98',Transit_Amplifying='#B84B46',Colonocyte='#7560A8',Goblet='#B58B20')
pal_donor<-c(DNW14='#2B6E98',DNW15='#B84B46',DNW9='#2E8079')
ink<-'#26343D';grey<-'#9AA4AD';blue<-'#356B91';red<-'#B84B46'
theme_set(theme_classic(base_size=7,base_family='sans')+theme(text=element_text(colour=ink),axis.text=element_text(colour=ink,size=7),axis.title=element_text(size=7),axis.line=element_line(linewidth=.18),axis.ticks=element_line(linewidth=.18),plot.title=element_text(size=8.5,face='bold'),plot.subtitle=element_text(size=7),legend.title=element_text(size=7),legend.text=element_text(size=6.7),legend.key.size=unit(2.5,'mm'),legend.margin=margin(0,0,0,0),strip.background=element_blank(),strip.text=element_text(size=7,face='bold'),plot.margin=margin(5,8,5,10),legend.spacing.x=unit(2,'mm'),legend.key.width=unit(4,'mm'),panel.spacing.x=unit(5,'mm')))
wrap<-function(x,n=25)vapply(x,function(z)paste(strwrap(z,n),collapse='\n'),character(1))
fmt<-function(x)ifelse(is.na(x),'NE',trimws(formatC(x,digits=3,format='g')))
modids<-c('CURATED_IFN_VISIBILITY','ICI_COLITIS_EPITHELIAL_UP_STRICT','ICI_COLITIS_EPITHELIAL_UP_NON_IFN','ICI_COLITIS_EPITHELIAL_LOSS_STRICT')
modlabs<-c('IFN reference','Disease-up','IFN-gene-excluded\nresidual','Disease-depleted')
modlab<-function(x)unname(setNames(modlabs,modids)[x])
families<-c('Absorptive epithelial','Immature epithelial','Mature absorptive epithelial','Secretory epithelial','Platform samples')
flabs<-c('Absorptive','Immature','Mature absorptive','Secretory','Platform samples')
ctlab<-c(Stem='Stem',Transit_Amplifying='TA',Colonocyte='Colonocyte',Goblet='Goblet')
pathlab<-function(x){z<-gsub('_',' ',sub('HALLMARK_','',x));z<-tools::toTitleCase(tolower(z));z<-gsub('Interferon Alpha Response','IFN-alpha response',z);z<-gsub('Interferon Gamma Response','IFN-gamma response',z);z<-gsub('Tnfa Signaling Via Nfkb','TNF / NF-kB',z,ignore.case=TRUE);z<-gsub('Epithelial Mesenchymal Transition','Epithelial-mesenchymal transition',z);z<-gsub('Il6 Jak Stat3 Signaling','IL6 / JAK / STAT3',z);z<-gsub('Il2 Stat5 Signaling','IL2 / STAT5',z);z<-gsub('G2m','G2M',z);z<-gsub('Myc','MYC',z);z<-gsub('Dna','DNA',z);z<-gsub('Uv','UV',z);z<-gsub('Mtorc1','mTORC1',z);z}
nes_scale<-function(d,name='NES') {lim<-max(abs(d),na.rm=TRUE);scale_fill_gradient2(low=blue,mid='white',high=red,midpoint=0,limits=c(-lim,lim),na.value='#D4D9DD',name=name)}
compose<-function(...,labs=NULL,ncol=1,rh=1,rw=1)plot_grid(...,labels=labs,ncol=ncol,rel_heights=rh,rel_widths=rw,label_size=9.5,label_fontface='bold')
