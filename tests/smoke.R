source('R/init.R')
stopifnot(!'DST'%in%paper_scope('G2','A'),! 'DST'%in%paper_scope('G2','B'),! 'DST'%in%paper_scope('G2','C'))
stopifnot(inherits(try(paper_check_methods('G2','DST','A'),silent=TRUE),'try-error'))
stopifnot(inherits(try(paper_check_methods('clinical','PR-GAM'),silent=TRUE),'try-error'))
stopifnot(length(paper_scope('clinical'))==17)
set.seed(14);n<-140;d<-data.frame(eid=paste0('s',1:n),row_id=1:n,Y=rexp(n,.04),delta=rbinom(n,1,.7),
 X1=rbinom(n,1,.4),X2=rbinom(n,1,.5),X3=rnorm(n),X4=rnorm(n),X5=rnorm(n),X6=rnorm(n),X7=rnorm(n))
s<-paper_gam(d,12,10,TRUE,TRUE);ctx<-paper_context(s,d,12,TRUE)
for(f in c('SPT','PR-GAM','Fixed-center','Uncalibrated','DST')) {
 z<-paper_candidate(ctx,.5,f);pred<-paper_predict(z,ctx,d)
 stopifnot(all(is.finite(pred)),all(pred>=.001 & pred<=.999),z$rho>=0,z$rho<=1,z$kkt<2e-6)
 if(f=='PR-GAM')stopifnot(z$rho==0)
 if(f=='Fixed-center')stopifnot(z$rho==1)
 if(f=='DST')stopifnot(z$rho==z$omega)
}
L<-cbind(rep(1,n),rep(1,n),rep(1,n),rep(1,n));stopifnot(paper_select(L)$tau==.5)
L[,1]<-.1;stopifnot(paper_select(L)$tau==0)
# The new core has no callable full-residual-channel estimator.
stopifnot(inherits(try(se_solver(ctx$B,ctx$c,ctx$r,ctx$y,ctx$w,ctx$q,.5,ctx$S,coef(ctx$anchor$fit),fixed_omega=NULL),silent=TRUE),'try-error'))
cat('Scope, formulas, displayed estimators and sole tau rule: PASS\n')
