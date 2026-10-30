source('R/init.R')
a<-commandArgs(TRUE);if(length(a)!=2)stop('Usage: Rscript prepare_realdata_sources.R CLEANED_DATA_DIR PRIVATE_SOURCE_CACHE')
d<-paper_clinical_data(a[1]);s<-d$source;dir.create(a[2],recursive=TRUE,showWarnings=FALSE)
key<-digest::digest(s,algo='sha256')
for(kind in c('GAM','Cox','NN'))for(h in if(kind=='Cox')24 else c(24,30,36)) {
 path<-file.path(a[2],paste0(kind,'_',if(kind=='Cox')'all'else h,'.rds'))
 if(file.exists(path)){stopifnot(readRDS(path)$data_sha256==key);next}
 fit<-switch(kind,GAM=paper_gam(s,h,k=10,source=TRUE,clinical=TRUE),Cox=paper_source_cox(s),NN=paper_nn(s,h,a[2]))
 saveRDS(list(data_sha256=key,fit=fit),path)
}
