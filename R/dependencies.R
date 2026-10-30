# Optional author implementations are loaded only for a requested displayed method.
load_transcox_sources<-function() {if(!exists('TransCox_R',mode='function'))source('R/transcox.R');invisible(TRUE)}
.score_load<-function(kind) {
 if(kind=='coxkl')stopifnot(requireNamespace('survkl',quietly=TRUE))
 if(kind=='potl' && !exists('.potl_api',inherits=TRUE)) {
  .potl_api<<-new.env(parent=globalenv())
  Rcpp::sourceCpp('external/POTL/EM.cpp',env=.potl_api,cacheDir=file.path(tempdir(),'spt_potl'),rebuild=FALSE)
 }
}
.nc_load<-function(kind) {
 if(kind=='discretekl')stopifnot(requireNamespace('DiscreteKL',quietly=TRUE)) else if(!exists('.aft_api',inherits=TRUE)) {
  suppressPackageStartupMessages(library(glmnet))
  .aft_api<<-new.env(parent=globalenv());sys.source('external/Trans-AFT/Simulation-Functions.R',.aft_api)
  .aft_api$cv.glmnet<-function(...) {
   args<-list(...);args$foldid<-sample(rep(seq_len(args$nfolds),length.out=nrow(args[[1]])))
   fit<-do.call(glmnet::cv.glmnet,args);.aft_api$last_cv<-fit;.aft_api$last_foldid<-args$foldid;fit
  }
 }
}
paper_load_coxtl<-function() {
 if(exists('.suite_coxtl_env',inherits=TRUE))return(invisible(NULL))
 .suite_coxtl_env<<-new.env(parent=globalenv())
 sys.source('external/CoxTL/run_CoxTL.R',.suite_coxtl_env);sys.source('external/CoxTL/density_opt.R',.suite_coxtl_env)
 .suite_coxtl_env$rcorr.cens<-Hmisc::rcorr.cens;.suite_coxtl_env$Surv<-survival::Surv;.suite_coxtl_env$strat<-rms::strat
 f<-.suite_coxtl_env$run_CoxTL;b<-body(f);stopifnot(identical(b[[length(b)]],quote(return(cox_final))))
 b[[length(b)]]<-quote({attr(cox_final,'native_tuning')<-list(scores=cv_record,folds=split_set,grid=lam_set,selected=lam_best,size_ratio=size_ratio);return(cox_final)})
 body(f)<-b;.suite_coxtl_env$run_CoxTL<-f
}
paper_log<-function(out,prefix) {
 e<-new.env();e$stage<-'init';e$write<-function(x,name)write.csv(x,file.path(out,paste0(prefix,'_',name,'.csv')),row.names=FALSE);e
}
paper_comparator<-function(method,s,t,te,hs,seed,out,clinical=FALSE) {
 oldopt<-options(spt.paper.clinical=clinical);on.exit(options(oldopt),add=TRUE)
 if(method=='CoxTL')paper_load_coxtl()
 if(method %in% c('Target-Cox','CoxKL','POTL-Cox')) {
  cfg<-jsonlite::fromJSON('config/score_transfer.json');cfg$horizons<-if(clinical)c(24,36,48)else hs
  if(method=='Target-Cox' && clinical) {
   f<-survival::coxph(reformulate(paste0('X',1:.ext_p),'survival::Surv(Y,delta)'),data=t,x=TRUE,y=TRUE)
   sf<-survival::survfit(f,newdata=te,se.fit=FALSE)
   return(clip(t(1-summary(sf,times=hs,extend=TRUE)$surv)))
  }
  if(method=='Target-Cox')o<-.score_source(t) else {
   sc<-.score_source(s)
   if(method=='CoxKL') {.score_load('coxkl');o<-if(clinical).score_coxkl(t,sc,seed,cfg)else .usim_coxkl_extended(t,sc,seed,cfg)}
   else {.score_load('potl');o<-.score_potl(t,sc,seed,cfg)}
   write.csv(o$trace,file.path(out,paste0(method,'_tuning.csv')),row.names=FALSE)
  }
  return(clip(.score_predict(o,.score_x(te),hs)))
 }
 if(clinical && method %in% c('DiscreteKL','Trans-AFT')) {
  cfg<-jsonlite::fromJSON('config/noncox.json');cfg$breaks<-seq(0,48,6)
  if(method=='DiscreteKL') {.nc_load('discretekl');sc<-.nc_disc_source(s,cfg);o<-.nc_discretekl(t,sc,seed,cfg)}
  else {.nc_load('transaft');o<-.nc_transaft(t,s,seed,cfg)}
  write.csv(o$trace,file.path(out,paste0(method,'_tuning.csv')),row.names=FALSE)
  return(clip(.nc_predict(o,te[paste0('X',1:.ext_p)],hs)))
 }
 .r100_extra_fit(method,s,t,te,hs,seed,paper_log(out,method))
}
