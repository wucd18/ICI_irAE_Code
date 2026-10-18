args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
O <- normalizePath(args[1],winslash='/')
hc3 <- function(y,X) {
 stopifnot(length(y)==nrow(X),all(is.finite(y)),all(is.finite(X)),qr(X)$rank==ncol(X))
 inv <- solve(crossprod(X));beta <- as.vector(inv%*%crossprod(X,y));res <- y-X%*%beta
 h <- rowSums((X%*%inv)*X);stopifnot(all(h<1-1e-8))
 v <- inv%*%crossprod(X,as.vector(res/(1-h))^2*X)%*%inv
 list(beta=beta,se=sqrt(diag(v)),res=as.vector(res),h=h,df=nrow(X)-ncol(X))
}
if (identical(Sys.getenv('ICI_FIXTURE_ONLY'),'1')) {
 x <- cbind(1,c(0,0,0,0,1,1,1,1),c(1,4,2,3,2,6,4,3))
 y <- c(2,1,4,3,7,6,8,9)
 a <- hc3(y,x);b <- lm.fit(x,y)
 stopifnot(max(abs(a$beta-b$coefficients))<1e-10,all(a$se>0),max(abs(hc3(2*y,x)$se-2*a$se))<1e-10)
 bad <- try(hc3(y,cbind(x,x[,2])),silent=TRUE);stopifnot(inherits(bad,'try-error'))
 cat('FIXTURE_PASS: coefficient reference, scaling, singular design rejection\n');quit(status=0)
}
stopifnot(file.exists(file.path(O,'review/conditional_analysis_plan_locked.json')))
d <- read.delim(file.path(O,'source_data/conditional_model_input.tsv'),check.names=FALSE)
outcomes <- c('ICI_COLITIS_EPITHELIAL_UP_NON_IFN','ICI_COLITIS_EPITHELIAL_LOSS_STRICT')
rows <- list();loo <- list();points <- list()
for (ds in unique(d$dataset)) for(ct in unique(d$celltype[d$dataset==ds])) {
 b <- d[d$dataset==ds & d$celltype==ct,];case <- if(ds=='GSE206300') 'irColitis' else 'CPI_colitis'
 disease <- as.integer(b$condition==case);ifn <- b$CURATED_IFN_VISIBILITY
 stopifnot(sum(disease==1)>=3,sum(disease==0)>=3,!anyDuplicated(b$patient_id))
 X <- cbind(intercept=1,disease=disease,IFN=ifn-mean(ifn));r2x <- summary(lm(disease~ifn))$r.squared
 for (outcome in outcomes) {
  y <- b[[outcome]];a <- hc3(y,X);u <- hc3(y,X[,1:2]);tval <- a$beta[2]/a$se[2];critical <- qt(.975,a$df)
  f0 <- lm(y~ifn);f1 <- lm(y~disease+ifn)
  rows[[length(rows)+1]] <- data.frame(dataset=ds,celltype=ct,module=outcome,n_case=sum(disease==1),n_control=sum(disease==0),unadjusted_beta=u$beta[2],unadjusted_se_HC3=u$se[2],adjusted_beta=a$beta[2],adjusted_se_HC3=a$se[2],ci_low=a$beta[2]-critical*a$se[2],ci_high=a$beta[2]+critical*a$se[2],p=2*pt(-abs(tval),a$df),residual_df=a$df,disease_IFN_VIF=1/(1-r2x),max_leverage=max(a$h),partial_R2=(deviance(f0)-deviance(f1))/deviance(f0),IFN_range_overlap=max(0,min(max(ifn[disease==1]),max(ifn[disease==0]))-max(min(ifn[disease==1]),min(ifn[disease==0]))),role=if(ds=='GSE206300')'independent_primary' else 'discovery_self_check')
  for(i in seq_len(nrow(b))) {
   aa <- hc3(y[-i],X[-i,,drop=FALSE])
   loo[[length(loo)+1]] <- data.frame(dataset=ds,celltype=ct,module=outcome,omitted_patient=b$patient_id[i],adjusted_beta=aa$beta[2])
  }
  points[[length(points)+1]] <- data.frame(dataset=ds,celltype=ct,module=outcome,patient_id=b$patient_id,condition=b$condition,IFN_score=ifn,oriented_score=y,IFN_adjusted_residual=residuals(f0),model_residual=a$res,leverage=a$h)
 }
}
r <- do.call(rbind,rows);l <- do.call(rbind,loo);pt <- do.call(rbind,points)
r$fdr <- NA_real_
for(ds in unique(r$dataset))r$fdr[r$dataset==ds] <- p.adjust(r$p[r$dataset==ds],method='BH')
for(i in seq_len(nrow(r))) {
 z <- l$adjusted_beta[l$dataset==r$dataset[i]&l$celltype==r$celltype[i]&l$module==r$module[i]]
 r$LOPO_beta_min[i] <- min(z);r$LOPO_beta_max[i] <- max(z);r$LOPO_all_expected_direction[i] <- all(z>0)
}
write.table(r,file.path(O,'source_data/conditional_results.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
write.table(l,file.path(O,'source_data/conditional_leave_patient_out.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
write.table(pt,file.path(O,'source_data/conditional_patient_diagnostics.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
capture.output(sessionInfo(),file=file.path(O,'review/conditional_R_sessionInfo.txt'))
print(r);cat('CONDITIONAL_ANALYSIS_COMPLETED\n')
