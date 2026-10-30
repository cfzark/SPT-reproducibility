# The paper estimator and ONLY its displayed matched ablations.
clip <- function(p) pmin(pmax(p,.001),.999)
.paper_cont <- function(clinical=FALSE) if(clinical)paste0('X',3:7) else paste0('X',3:.ext_p)
paper_formula <- function(k=5,interaction=FALSE,source=FALSE,clinical=FALSE,group=FALSE,response='D_t') {
 if(!clinical) f <- .gcv_formula(interaction,k,source,group=group) else {
  rhs<-c('X2','X8',sprintf("s(X%d_sc,k=%d,bs='ts')",3:7,k),if(group)'target_group')
  f<-as.formula(paste('D_t~',paste(rhs,collapse='+')))
 }
 f[[2]]<-as.name(response);f
}
paper_gam <- function(d,h,k=5,interaction=FALSE,source=FALSE,clinical=FALSE,choose=FALSE) {
 aux<-make_Dt_weights(d,h);sc<-make_scaler(d,.paper_cont(clinical));z<-apply_scaler(d,sc)
 z$D_t<-aux$D_t;z$w_ipcw<-aux$weights
 if(choose)interaction<-.gcv_choose(z,z$w_ipcw,k,if(source)'source'else'target',h,source)$interaction
 f<-paper_formula(k,interaction,source,clinical)
 o<-mgcv::gam(f,data=z,family=quasibinomial(),weights=w_ipcw,method='REML',gamma=1.4,select=FALSE)
 o<-.cur_retry_gam(o);stopifnot(isTRUE(o$converged),all(is.finite(coef(o))))
 list(fit=o,scaler=sc,kind='gam',interaction=interaction,clinical=clinical,k=k)
}
paper_source_lp <- function(o,d,h=NULL) {
 if(o$kind=='nn') return(qlogis(clip(.inn_prob(o,d))))
 z<-apply_scaler(d,o$scaler)
 if(o$kind=='cox') {
  lp<-as.numeric(predict(o$fit,z,type='lp',reference='zero'))
  H<-c(0,o$baseline$hazard)[findInterval(h,o$baseline$time)+1L]
  return(qlogis(clip(-expm1(-H*exp(lp)))))
 }
 as.numeric(predict(o$fit,z,type='link'))
}
paper_context <- function(src,d,h,interaction=FALSE,choose=FALSE,clinical=FALSE) {
 a<-paper_gam(d,h,interaction=interaction,choose=choose,clinical=clinical)
 z<-apply_scaler(d,a$scaler);z$source_lp<-paper_source_lp(src,d,h)
 pr<-mgcv::gam(paper_formula(5,a$interaction,FALSE,clinical,response='source_lp'),data=z,
  family=gaussian(),method='REML',gamma=1.4)
 B<-unname(predict(a$fit,z,type='lpmatrix'));PB<-unname(predict(pr,z,type='lpmatrix'))
 stopifnot(isTRUE(pr$converged),identical(dim(B),dim(PB)),max(abs(B-PB))<1e-10,qr(B,tol=1e-10)$rank==ncol(B))
 c<-unname(coef(pr));aux<-make_Dt_weights(d,h)
 list(anchor=a,projection=pr,source_obj=src,horizon=h,B=B,c=c,r=z$source_lp-drop(B%*%c),
      S=paper_penalty(a$fit),y=aux$D_t,w=aux$weights,q=clip(plogis(z$source_lp)),source_lp=z$source_lp)
}
paper_candidate <- function(ctx,tau,family) {
 stopifnot(family %in% c('SPT','PR-GAM','Fixed-center','Uncalibrated','DST'))
 cal<-uro_calibrate_teacher(ctx$source_lp,ctx$y,ctx$w)
 if(family=='DST')z<-usi_single_fit(ctx,cal,tau) else {
  q<-if(family=='Uncalibrated'||tau==0)ctx$q else cal$q
  rho<-switch(family,'PR-GAM'=0,'Fixed-center'=1,NULL)
  z<-se_solver(ctx$B,ctx$c,ctx$r,ctx$y,ctx$w,q,tau,ctx$S,coef(ctx$anchor$fit),fixed_rho=rho,fixed_omega=0)
 }
 z$tau<-tau;z$family<-family;z
}
paper_predict <- function(fit,ctx,d) {
 B<-unname(predict(ctx$anchor$fit,apply_scaler(d,ctx$anchor$scaler),type='lpmatrix'))
 eta<-drop(B%*%fit$theta)
 if(fit$family=='DST') eta<-eta+fit$omega*(paper_source_lp(ctx$source_obj,d,ctx$horizon)-drop(B%*%ctx$c))
 clip(plogis(eta))
}
paper_select <- function(loss,taus=c(0,.1,.5,1)) {
 risk<-colMeans(loss);j<-order(risk,taus)[1];diff<-sweep(loss,1,loss[,j],'-')
 se<-apply(diff,2,sd)/sqrt(nrow(loss));tol<-64*.Machine$double.eps*max(1,max(abs(risk)))
 eligible<-colMeans(diff)<=se+tol;ref<-which(taus==.5);stopifnot(length(ref)==1,all(is.finite(loss)))
 selected<-if(eligible[ref])ref else j
 list(index=selected,tau=taus[selected],table=data.frame(tau=taus,risk=risk,difference=colMeans(diff),paired_se=se,eligible=eligible))
}
paper_transfer <- function(src,tr,te,h,folds,labels,interaction=FALSE,choose=FALSE,clinical=FALSE) {
 # labels map to family + fixed tau; all other labels use the sole published selector.
 stopifnot(all(labels %in% c('SPT','PR-GAM','Fixed-center','Uncalibrated','DST','SPT-tau0','DST-tau0')))
 family<-sub('-tau0$','',labels);auto<-!grepl('-tau0$',labels);taus<-c(0,.1,.5,1)
 losses<-setNames(lapply(unique(family[auto]),function(x)matrix(NA_real_,nrow(tr),4)),unique(family[auto]))
 if(any(auto))for(k in sort(unique(folds))) {
  ct<-tr[folds!=k,];va<-tr[folds==k,];ctx<-paper_context(src,ct,h,interaction,choose,clinical)
  km<-.guard_km_fit(ct)
  for(f in names(losses)) {
   pp<-sapply(taus,function(t)paper_predict(paper_candidate(ctx,t,f),ctx,va))
   losses[[f]][folds==k,]<-ucv_validation(va,h,pp,km)$loss
  }
 }
 ctx<-paper_context(src,tr,h,interaction,choose,clinical);fits<-selections<-list()
 for(i in seq_along(labels)) {
  sel<-if(auto[i])paper_select(losses[[family[i]]]) else list(tau=0)
  fits[[labels[i]]]<-paper_candidate(ctx,sel$tau,family[i]);selections[[labels[i]]]<-sel
 }
 pp<-sapply(fits,paper_predict,ctx=ctx,d=te)
 pp<-cbind('T-GAM'=clip(plogis(paper_source_lp(ctx$anchor,te))),pp)
 list(predictions=pp,context=ctx,fits=fits,selections=selections)
}
paper_pool <- function(s,t,te,h,group=FALSE,clinical=FALSE) {
 sa<-make_Dt_weights(s,h);ta<-make_Dt_weights(t,h);d<-rbind(s,t)
 sc<-make_scaler(d,.paper_cont(clinical));d<-apply_scaler(d,sc)
 d$D_t<-c(sa$D_t,ta$D_t);d$w<-c(sa$weights,ta$weights)
 d$target_group<-c(rep(0,nrow(s)),rep(1,nrow(t)))
 f<-paper_formula(5,!clinical,TRUE,clinical,group)
 o<-mgcv::gam(f,data=d,family=quasibinomial(),weights=w,method='REML',gamma=1.4,select=FALSE)
 stopifnot(isTRUE(o$converged),all(is.finite(coef(o))))
 z<-apply_scaler(te,sc);z$target_group<-1;clip(as.numeric(predict(o,z,type='response')))
}
paper_metric <- function(d,h,p,q=NULL) {
 p<-clip(p);y<-as.integer(d$Y<=h & d$delta==1);w<-.frozen_km_weights(d$Y,d$delta,h)
 out<-data.frame(brier=mean(w*(y-p)^2),auc=.usim_auc(p,w*y,w*(d$Y>h)),
 uno=survival::concordance(survival::Surv(Y,delta)~p,data=d,reverse=TRUE,timewt='n/G2',ymax=h)$concordance)
 if(!is.null(q)) {
  out$true_brier<-mean((p-q)^2+q*(1-q));out$true_auc<-.usim_auc(p,q,1-q)
  out$kl_risk<-mean(-q*log(p)-(1-q)*log1p(-p)-.usim_entropy(q))
 }
 out
}
