# Seven-variable UKB B -> SPRINT standard; death, 24/30/36 months only.
.clean_map<-c(X2='SM_history',X3='age',X4='bmi',X5='GFR',X6='cholesterol',X7='creatinine',X8='sex')
source('R/nn.R')
paper_clinical_data<-function(root) {
 read<-function(name,prefix,n) {
  a<-read.csv(file.path(root,paste0(name,'.csv')),colClasses=c(eid='character'))
  stopifnot(nrow(a)==n,all(c('eid','follow_up_time','death_event',unname(.clean_map))%in%names(a)),!anyDuplicated(a$eid))
  d<-data.frame(eid=paste(prefix,a$eid,sep=':'),Y=a$follow_up_time,delta=a$death_event,row_id=seq_len(nrow(a)))
  for(v in names(.clean_map))d[[v]]<-a[[.clean_map[[v]]]]
  stopifnot(all(is.finite(as.matrix(d[c('Y','delta',names(.clean_map))]))),all(d$Y>0),all(d$delta%in%0:1))
  if('followup_days'%in%names(a))stopifnot(max(abs(d$Y-a$followup_days/30.4375))<1e-10)
  d
 }
 list(source=read('ukb_B_primary','UKB',7059),target=read('sprint_standard','SPRINT',1126))
}
paper_nn<-function(d,h,out) {
 ff<-.paper_folds(d,36,5,2026092601L);grid<-expand.grid(size=c(4L,8L),decay=c(.1,1,10));L<-matrix(NA_real_,nrow(d),nrow(grid))
 for(g in seq_len(nrow(grid)))for(k in 1:5) {
  tr<-d[ff!=k,];va<-d[ff==k,];o<-.inn_fit(tr,h,grid$size[g],grid$decay[g],20261012L+h*1000L+k*10L+1:3)
  L[ff==k,g]<-ucv_validation(va,h,matrix(.inn_prob(o,va),ncol=1),.guard_km_fit(tr))$loss
 }
 grid$loss<-colMeans(L);g<-which.min(grid$loss);grid$selected<-seq_len(nrow(grid))==g
 write.csv(grid,file.path(out,paste0('NN_grid_',h,'.csv')),row.names=FALSE)
 o<-.inn_fit(d,h,grid$size[g],grid$decay[g],20261012L+h*1000L+1:3);o$kind<-'nn';o
}
paper_source_cox<-function(s) {
 sc<-make_scaler(s,paste0('X',3:7));z<-apply_scaler(s,sc)
 f<-reformulate(c('X2','X8',paste0('X',3:7,'_sc')),'survival::Surv(Y,delta)')
 o<-survival::coxph(f,data=z,ties='efron',x=TRUE,y=TRUE,model=TRUE,singular.ok=FALSE,control=survival::coxph.control(iter.max=100,eps=1e-9))
 stopifnot(all(is.finite(coef(o))),max(o$iter)<100)
 list(fit=o,scaler=sc,baseline=survival::basehaz(o,centered=FALSE),kind='cox')
}
paper_recalibrate<-function(lp,y,w,newlp) {
 mu<-sum(w*lp)/sum(w);ss<-sqrt(sum(w*(lp-mu)^2)/sum(w));stopifnot(ss>0)
 x<-cbind(1,(lp-mu)/ss);o<-glm.fit(x,y,weights=w,family=quasibinomial(),control=glm.control(maxit=200,epsilon=1e-10))
 stopifnot(isTRUE(o$converged),all(is.finite(o$coefficients)),max(abs(drop(crossprod(x,w*(o$fitted.values-y)))))/sum(w)<1e-6)
 clip(plogis(drop(cbind(1,(newlp-mu)/ss)%*%o$coefficients)))
}
paper_relabel<-function(d) {x<-d[names(.clean_map)];names(x)<-paste0('X',1:7);cbind(d[c('eid','row_id','Y','delta')],x)}
paper_clinical<-function(data_root,rep,out,source_cache,methods=NULL,retain=FALSE) {
 stopifnot(rep%in%1:100);.ext_p<<-7L
 if(is.null(methods))methods<-paper_scope('clinical');paper_check_methods('clinical',methods)
 dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(source_cache,recursive=TRUE,showWarnings=FALSE)
 if(file.exists(file.path(out,'metrics.csv')))stop('Output already exists')
 d<-paper_clinical_data(data_root);s<-d$source;t<-d$target;hs<-c(24,30,36)
 ff<-.paper_folds(t,36,5,20260816+10007L*rep)
 write.csv(data.frame(row_id=t$row_id,fold=ff),file.path(out,'outer_folds.csv'),row.names=FALSE)
 getsource<-function(kind,h) {
  path<-file.path(source_cache,paste0(kind,'_',if(kind=='Cox')'all'else h,'.rds'))
  key<-digest::digest(s,algo='sha256')
  if(file.exists(path)){z<-readRDS(path);stopifnot(identical(z$data_sha256,key));return(z$fit)}
  fit<-switch(kind,GAM=paper_gam(s,h,k=10,source=TRUE,clinical=TRUE),Cox=paper_source_cox(s),NN=paper_nn(s,h,out))
  saveRDS(list(data_sha256=key,fit=fit),path);fit
 }
 P<-array(NA_real_,c(nrow(t),3,length(methods)),dimnames=list(NULL,NULL,methods));dec<-list()
 core<-intersect(methods,c('SPT-GAM','SPT-Cox','SPT-NN'))
 for(j in 1:3) {
  h<-hs[j]
  if(any(c('Source-GAM','SPT-GAM','Score-recalibration')%in%methods))gs<-getsource('GAM',h)
  if('Source-GAM'%in%methods)P[,j,'Source-GAM']<-clip(plogis(paper_source_lp(gs,t,h)))
  for(f in 1:5) {
   tr<-t[ff!=f,];te<-t[ff==f,];inner<-.paper_folds(tr,36,3,20260922+10007L*rep+101L*f)
   for(m in core) {
    src<-getsource(sub('SPT-','',m),h);z<-paper_transfer(src,tr,te,h,inner,'SPT',clinical=TRUE)
    P[ff==f,j,m]<-z$predictions[,'SPT']
    if('T-GAM'%in%methods)P[ff==f,j,'T-GAM']<-z$predictions[,'T-GAM']
    dec[[length(dec)+1L]]<-data.frame(method=m,horizon=h,fold=f,tau=z$fits$SPT$tau,rho=z$fits$SPT$rho)
    if(retain)saveRDS(z,file.path(out,paste0('private_',m,'_h',h,'_fold',f,'.rds')))
   }
   if(!length(core)&&'T-GAM'%in%methods)P[ff==f,j,'T-GAM']<-clip(plogis(paper_source_lp(paper_gam(tr,h,clinical=TRUE),te)))
   for(m in intersect(methods,c('Pooled-GAM','Stratified-GAM')))P[ff==f,j,m]<-paper_pool(s,tr,te,h,m=='Stratified-GAM',TRUE)
   if('Score-recalibration'%in%methods){a<-make_Dt_weights(tr,h);P[ff==f,j,'Score-recalibration']<-paper_recalibrate(paper_source_lp(gs,tr,h),a$D_t,a$weights,paper_source_lp(gs,te,h))}
  }
 }
 extra<-setdiff(methods,c(core,'T-GAM','Source-GAM','Pooled-GAM','Stratified-GAM','Score-recalibration'))
 for(f in 1:5)for(m in extra) {
  base<-if(m%in%c('CoxKL','POTL-Cox'))20260916L else if(m%in%c('DiscreteKL','Trans-AFT'))20260917L else 20260915L
  subout<-file.path(out,paste0('fold',f,'_',m));dir.create(subout,showWarnings=FALSE)
  P[ff==f,,m]<-paper_comparator(m,paper_relabel(s),paper_relabel(t[ff!=f,]),paper_relabel(t[ff==f,]),hs,base+10007L*rep+101L*f,subout,TRUE)
 }
 stopifnot(all(is.finite(P)))
 met<-do.call(rbind,lapply(methods,function(m)do.call(rbind,lapply(1:3,function(j){
  p<-P[,j,m];mm<-paper_metric(t,hs[j],p)
  mm$auc<-tail(timeROC::timeROC(T=t$Y,delta=t$delta,marker=qlogis(p),cause=1,weighting='marginal',times=hs[j],iid=FALSE)$AUC,1)
  mm$brier<-mean(estimate_weights(t,hs[j],w_cap=Inf)*(as.integer(t$Y<=hs[j]&t$delta==1)-p)^2)
  cbind(data.frame(design='clinical',replication=rep,horizon=hs[j],method=m),mm)
 }))))
 write.csv(met,file.path(out,'metrics.csv'),row.names=FALSE)
 if(length(dec))write.csv(do.call(rbind,dec),file.path(out,'selection.csv'),row.names=FALSE)
 if(retain)saveRDS(list(predictions=P,folds=ff,ids=t$eid,time=t$Y,event=t$delta),file.path(out,'private_predictions.rds'))
 capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'));invisible(met)
}
