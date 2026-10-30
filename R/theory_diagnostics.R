inv<-function(H)chol2inv(chol(H))
km_if<-function(t,event,h,f,p,w) {
 n<-length(t);a<-w*f*(as.numeric(t<=h & event==1)-p)
 times<-sort(unique(t));idx<-match(t,times);ng<-length(times)
 nn<-tabulate(idx,nbins=ng);dc<-tabulate(idx[event==0],nbins=ng)
 ag<-as.numeric(rowsum(a,idx,reorder=TRUE));stopifnot(length(ag)==ng)
 R<-rev(cumsum(rev(nn)));ta<-rev(cumsum(rev(ag)))
 eligible<-dc>0 & times<=h;b<-numeric(ng)
 stopifnot(all(R[eligible]>dc[eligible]));b[eligible]<-ta[eligible]/(R[eligible]-dc[eligible])
 corr<-as.numeric(event==0 & t<=h)*b[idx]-cumsum(b*dc/R)[idx]
 psi<-a+corr;list(psi=psi-mean(psi),nu=var(psi),nu_fixed=var(a),mean=mean(a))
}
diag_one <- function(ctx,Bte,rte,theta,wt,p_tr,rho,t,truthq) {
 B<-ctx$B;n<-nrow(B);A<-ctx$S/sqrt(n);c<-ctx$c
 H<-crossprod(Bte,Bte*wt)/nrow(Bte);Hi<-inv(H)
 h<-drop(crossprod(Bte,wt*rte)/nrow(Bte));alpha<-drop(Hi%*%h)
 s<-mean(wt*rte^2)-sum(h*alpha);stopifnot(s>0)
 f<-ctx$r-drop(B%*%alpha)
 ii<-km_if(t$Y,t$delta,ctx$horizon,f,p_tr,ctx$w)
 nu<-ii$nu;sd<-sqrt(nu);a<-drop(A%*%theta);ac<-drop(A%*%c)
 k<-sum((c+alpha)*a)
 GK<-rho*sum(ac*(Hi%*%a))-.5*rho^2*sum(ac*(Hi%*%ac))
 z<-k/sd;D<-((nu+k*k)*pnorm(z)+k*sd*dnorm(z))/(2*s)
 pzero<-pnorm(-z)
 Hhalf<-chol(H);L<-backsolve(Hhalf,diag(ncol(H)))
 rel<-crossprod(L,(ctx$S/n)%*%L)
 ev<-eigen(rel,symmetric=TRUE,only.values=TRUE)$values
 list(GK=GK,s=s,k=k,nu=nu,nu_fixed=ii$nu_fixed,D=D,total=GK+D,pzero=pzero,
      rho_formula=if(sum(c*ac)>0)max(0,min(1,sum(c*a)/sum(c*ac))) else 0,
      A_norm=max(eigen(A,symmetric=TRUE,only.values=TRUE)$values),
      penalty_curvature_max=max(ev),penalty_curvature_median=median(ev),
      score_eff=sqrt(n)*ii$mean,t_local=max(0,sqrt(n)*ii$mean+k)/s,
      Hi=Hi,H=H,alpha=alpha,if_score=ii$psi)
}
