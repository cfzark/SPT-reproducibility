paper_draw<-function(cell,rep) {
 .ext_p<<-cell$p;seed<-1012610042L+1000L*rep;setting<-match(cell$setting,LETTERS)
 s<-.sg_draw(3000,setting,'source',cell$distribution,cell$sigma,cell$source_lambda,seed+1,'source')
 t<-.sg_draw(cell$n_target,setting,'target',cell$distribution,cell$sigma,cell$target_lambda,seed+2,'target')
 te<-.sg_draw(5000,setting,'target',cell$distribution,cell$sigma,cell$test_lambda,seed+3,'test')
 hs<-as.numeric(quantile(t$data$Y,c(.25,.5,.75),type=7))
 te<-.sg_truth(te,hs,setting,'target',cell$distribution,cell$sigma)
 list(source=s,target=t,test=te,horizons=hs,seed=seed)
}
paper_stress_g<-function(x,dist,sigma,kind) {
 g<-.ext_shared_g(x,1,'source')
 switch(kind,level=g+sigma,shape=g+.8*(x$X4^2-1)/sqrt(.8)-.8*x$X5,reverse=-g)
}
paper_stress_source<-function(z,cell,kind) {
 # Calibrate only the source censoring distribution, using the historical fixed population seed.
 set.seed(810081L);x<-.usim_x(100000,'source');u<-runif(100000)
 tt<-exp(log(36)-paper_stress_g(x,cell$distribution,cell$sigma,kind)+cell$sigma*.usim_eps(u,cell$distribution))
 rate<-exp(uniroot(function(l)mean(-expm1(-exp(l)*tt))-.4,c(-50,30),tol=1e-11)$root)
 d<-z$source;oldg<-.ext_shared_g(d$data,1,'source');newg<-paper_stress_g(d$data,cell$distribution,cell$sigma,kind)
 d$event_time<-d$event_time*exp(oldg-newg);d$censor_time<-d$censor_time*cell$source_lambda/rate
 d$data$Y<-pmin(d$event_time,d$censor_time);d$data$delta<-as.integer(d$event_time<=d$censor_time)
 z$source<-d;z
}
paper_simulation<-function(design,cell_id,rep,out,methods=NULL,n_target=NULL,retain=FALSE) {
 stopifnot(rep %in% 1:100)
 cells<-read.csv('config/simulation_cells.csv');cell<-cells[cells$cell==cell_id,,drop=FALSE];stopifnot(nrow(cell)==1)
 if(design %in% c('G2','G3','G5'))stopifnot(cell$group==design)
 if(design=='structural_D')stopifnot(cell$group=='G2',cell$setting=='D')
 if(design %in% c('level','shape','reverse'))stopifnot(cell$group=='G2',cell$setting=='A')
 if(design=='large_n') {stopifnot(cell$setting=='D',cell$distribution=='loglogistic',n_target %in% c(500,1000,2000,4000,10000));cell$n_target<-n_target}
 if(is.null(methods))methods<-paper_scope(design,cell$setting)
 paper_check_methods(design,methods,cell$setting)
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 if(file.exists(file.path(out,'metrics.csv')))stop('Output already exists; use a separate directory')
 .gcv_trace<<-list();.cens_fallbacks<<-list()
 z<-paper_draw(cell,rep);if(design %in% c('level','shape','reverse'))z<-paper_stress_source(z,cell,design)
 s<-z$source$data;t<-z$target$data;te<-z$test$data;hs<-z$horizons
 choose<-design %in% c('G2','structural_D');interaction<-cell$setting!='D'
 inner<-.paper_folds(t,hs[2],3,z$seed+5L)
 write.csv(data.frame(row_id=t$row_id,fold=inner),file.path(out,'inner_folds.csv'),row.names=FALSE)
 pp<-array(NA_real_,c(nrow(te),3,length(methods)),dimnames=list(NULL,NULL,methods));fits<-list()
 core<-intersect(methods,c('SPT','PR-GAM','Fixed-center','Uncalibrated','DST','SPT-tau0','DST-tau0'))
 for(j in 1:3) {
  h<-hs[j];needsource<-length(core)||'Source-GAM'%in%methods
  if(needsource)src<-paper_gam(s,h,10,TRUE,TRUE,choose=choose)
  if(length(core)) {
   o<-paper_transfer(src,t,te,h,inner,core,interaction,choose)
   for(m in intersect(colnames(o$predictions),methods))pp[,j,m]<-o$predictions[,m]
   if(retain)fits[[as.character(j)]]<-o
   write.csv(do.call(rbind,lapply(names(o$fits),function(m)data.frame(method=m,tau=o$fits[[m]]$tau,rho=o$fits[[m]]$rho,omega=o$fits[[m]]$omega))),file.path(out,paste0('selection_Q',j,'.csv')),row.names=FALSE)
  } else if('T-GAM'%in%methods)pp[,j,'T-GAM']<-clip(plogis(paper_source_lp(paper_gam(t,h,interaction=interaction,choose=choose),te)))
  if('Source-GAM'%in%methods)pp[,j,'Source-GAM']<-clip(plogis(paper_source_lp(src,te)))
  for(m in intersect(methods,c('Pooled-GAM','Stratified-GAM')))pp[,j,m]<-paper_pool(s,t,te,h,m=='Stratified-GAM')
  if('Oracle'%in%methods)pp[,j,'Oracle']<-clip(z$test$true_risk[,j])
 }
 extra<-setdiff(methods,c(core,'T-GAM','Source-GAM','Pooled-GAM','Stratified-GAM','Oracle'))
 for(m in extra) {
  seed<-z$seed+if(m=='CoxKL')4L else if(m=='POTL-Cox')14L else 5L
  pp[,,m]<-paper_comparator(m,s,t,te,hs,seed,out)
 }
 stopifnot(all(is.finite(pp)))
 met<-do.call(rbind,lapply(methods,function(m)do.call(rbind,lapply(1:3,function(j)cbind(data.frame(design=design,setting=cell$setting,distribution=cell$distribution,n_target=cell$n_target,p=cell$p,replication=rep,quantile=c(.25,.5,.75)[j],horizon=hs[j],method=m),paper_metric(te,hs[j],pp[,j,m],z$test$true_risk[,j]))))))
 # KL is reported only by the large-n figure; DGP checks use Brier/AUC.
 if(design!='large_n')met$kl_risk<-NULL
 write.csv(met,file.path(out,'metrics.csv'),row.names=FALSE)
 if(length(.gcv_trace))write.csv(do.call(rbind,.gcv_trace),file.path(out,'structure.csv'),row.names=FALSE)
 if(length(.cens_fallbacks))write.csv(do.call(rbind,.cens_fallbacks),file.path(out,'numerical_fallbacks.csv'),row.names=FALSE)
 if(retain)saveRDS(list(data=z,fits=fits,predictions=pp),file.path(out,'simulation_models.rds'))
 capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'));invisible(met)
}
