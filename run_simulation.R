source('R/init.R')
a<-commandArgs(TRUE)
if(length(a)<4)stop('Usage: Rscript run_simulation.R DESIGN CELL REP OUTPUT [comma-separated-methods|all] [large_n_N]')
methods<-if(length(a)<5||a[5]=='all')NULL else strsplit(a[5],',',fixed=TRUE)[[1]]
paper_simulation(a[1],as.integer(a[2]),as.integer(a[3]),a[4],methods,
 n_target=if(length(a)>=6)as.integer(a[6])else NULL,retain=identical(Sys.getenv('SPT_SAVE_MODELS'),'1'))
