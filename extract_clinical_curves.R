source('R/init.R')
a<-commandArgs(TRUE);if(length(a)!=2)stop('Usage: Rscript extract_clinical_curves.R PRIVATE_REP_OUTPUT_ROOT OUTPUT.csv')
g<-read.csv('reference/clinical_curve_path.csv');ref<-read.csv('reference/clinical_curve_reference.csv')
x<-g[names(.clean_map)];r<-x[1,,drop=FALSE]
r$X2<-0;r$X3<-ref$age;r$X4<-ref$bmi;r$X5<-ref$egfr;r$X6<-ref$cholesterol;r$X7<-ref$creatinine;r$X8<-1
allx<-rbind(r,x);rows<-list()
for(rep in 1:100)for(fold in 1:5) {
 path<-file.path(a[1],sprintf('rep_%03d',rep),paste0('private_SPT-GAM_h36_fold',fold,'.rds'))
 z<-readRDS(path);ctx<-z$context;B<-predict(ctx$anchor$fit,apply_scaler(allx,ctx$anchor$scaler),type='lpmatrix')
 b<-sweep(B[-1,,drop=FALSE],2,B[1,],'-')
 v<-b%*%cbind('T-GAM'=coef(ctx$anchor$fit),'Source projection'=ctx$c,SPT=z$fits$SPT$theta)
 sl<-paper_source_lp(ctx$source_obj,allx,36);v<-cbind(v,'Source GAM'=sl[-1]-sl[1])
 rows[[length(rows)+1L]]<-cbind(data.frame(replication=rep,fold=fold),g[c('panel','x')],as.data.frame(v))
}
long<-data.table::melt(data.table::rbindlist(rows),id.vars=c('replication','fold','panel','x'),variable.name='method',value.name='value')
ans<-long[,list(mean=mean(value)),by=list(panel,x,method)];data.table::fwrite(ans,a[2])
