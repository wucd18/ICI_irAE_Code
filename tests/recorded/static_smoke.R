args<-commandArgs(TRUE)
p<-file.path(args[1],'presentation_code',c('plot_v1.R','plot_spatial_extension.R','plot_xenium_extension.R'))
expr<-as.expression(unlist(lapply(p,parse),recursive=FALSE))
calls<-character()
walk<-function(x){if(missing(x))return();if(is.call(x)){calls<<-c(calls,as.character(x[[1]])[1]);for(y in as.list(x)[-1])walk(y)}else if(is.expression(x)||is.pairlist(x))for(y in x)walk(y)}
walk(expr)
forbidden<-c('source','sys.source','system','system2','DESeq','glm','lm','coxph','rma','fgseaSimple','fgseaMultilevel','install.packages','download.file','unlink','file.remove')
stopifnot(!any(calls%in%forbidden))
fixture<-data.frame(NES=c(-2,1,NA),direction=c('down','up','down'),padj=c(.03,.4,NA))
oriented<-ifelse(fixture$direction=='down',-fixture$NES,fixture$NES)
stopifnot(identical(oriented,c(2,1,NA_real_)),is.na(oriented[3]))
# A saved null result and an absent result are distinct.
mark<-ifelse(is.na(fixture$NES),'NE',ifelse(fixture$padj<.05,'*',''))
stopifnot(identical(mark,c('*','','NE')))
# Python booleans must not become character row-name indexing in R.
flag<-c('True','False','TRUE','FALSE');take<-match(tolower(flag),c('false','true'))-1L
stopifnot(identical(which(take==1L),c(1L,3L)),!anyNA(take))
# Pure coordinate tiling retains every point, including both sides of a gap.
yy<-c(1,2,3,50,51);gap<-diff(yy);cut<-mean(yy[which.max(gap)+0:1]);tiles<-ifelse(yy<=cut,1L,2L)
stopifnot(identical(tiles,c(1L,1L,1L,2L,2L)),sum(table(tiles))==length(yy))
cat('PASS: R parse; prohibited-call audit; orientation and missingness fixtures.\n')
