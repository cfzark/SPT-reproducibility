source('R/init.R');source('R/theory_diagnostics.R')
a<-commandArgs(TRUE);if(length(a)!=2)stop('Usage: Rscript diagnose_large_n.R simulation_models.rds OUTPUT.csv')
z<-readRDS(a[1]);tr<-z$data$target$data;te<-z$data$test$data;n<-nrow(tr);rows<-list()
for(j in 1:3) {
 o<-z$fits[[as.character(j)]];ctx<-o$context
 stopifnot(!ctx$anchor$interaction,all(c('SPT-tau0','DST-tau0')%in%names(o$fits)))
 Bte<-unname(predict(ctx$anchor$fit,apply_scaler(te,ctx$anchor$scaler),type='lpmatrix'))
 rt<-paper_source_lp(ctx$source_obj,te,ctx$horizon)-drop(Bte%*%ctx$c)
 theta<-coef(ctx$anchor$fit);pt<-plogis(drop(ctx$B%*%theta));pe<-plogis(drop(Bte%*%theta));q<-z$data$test$true_risk[,j]
 dd<-diag_one(ctx,Bte,rt,theta,pe*(1-pe),pt,o$fits[['SPT-tau0']]$rho,tr,q)
 pp<-o$predictions;ce<-function(p)mean(-q*log(clip(p))-(1-q)*log1p(-clip(p)))
 rows[[j]]<-data.frame(n=n,horizon_index=j,A_norm=dd$A_norm,D=dd$D,k=dd$k,s=dd$s,nu=dd$nu,obs_DST_TG=n*(ce(pp[,'DST-tau0'])-ce(pp[,'T-GAM'])))
}
write.csv(do.call(rbind,rows),a[2],row.names=FALSE)
