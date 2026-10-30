source('R/init.R')
a<-commandArgs(TRUE)
if(length(a)<4)stop('Usage: Rscript run_realdata.R CLEANED_DATA_DIR REP OUTPUT PRIVATE_SOURCE_CACHE [comma-separated-methods|all]')
methods<-if(length(a)<5||a[5]=='all')NULL else strsplit(a[5],',',fixed=TRUE)[[1]]
paper_clinical(a[1],as.integer(a[2]),a[3],a[4],methods,retain=identical(Sys.getenv('SPT_SAVE_MODELS'),'1'))
