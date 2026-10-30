library(data.table)
out <- 'output';dir.create(out,showWarnings=FALSE)
z <- fread('reference/clinical_curves.csv')

cols<-c('T-GAM'='#595959','Source GAM'='#009E73','Source projection'='#CC79A7','SPT'='#0072B2')
methods<-names(cols);types<-c('T-GAM'=1,'Source GAM'=4,'Source projection'=3,'SPT'=1)
render<-function() {
 layout(matrix(c(1,2,3,4,5,5),3,2,byrow=TRUE),heights=c(1,1,.14))
 par(oma=c(0,2.4,0,0),mar=c(4.0,3.6,1.0,1.0),mgp=c(2.5,.65,0),tcl=-.25,family='sans',cex=.94)
 for(pname in c('Age','eGFR','Ever smoked','Sex')) {
  d<-z[z$panel==pname,];binary<-pname %in% c('Ever smoked','Sex')
  yl<-range(c(0,d$mean));pad<-max(diff(yl)*.08,.035);yl<-yl+c(-pad,pad)
  xr<-if(binary)c(-.08,1.08) else range(d$x)
  xl<-switch(pname,Age='Age (years)',eGFR=expression(eGFR~(mL/min/1.73~m^2)),pname)
  plot(NA,xlim=xr,ylim=yl,xlab=xl,ylab='',axes=FALSE)
  axis(2,at=pretty(yl,n=if(pname=='Sex')3 else 4),las=1);if(binary)axis(1,at=c(0,1),labels=if(pname=='Sex')c('Female','Male') else c('No','Yes')) else axis(1)
  box(col='#B0B0B0');abline(h=0,col='#D0D0D0',lwd=.8)
  for(m in methods) {
   q<-d[as.character(d$method)==m,];q<-q[order(x)]
   lines(q$x,q$mean,col=cols[m],lty=types[m],lwd=2)
   if(binary)points(q$x,q$mean,col=cols[m],pch=16,cex=.8)
  }
 }
 par(mar=c(0,0,0,0));plot.new()
 legend('center',legend=methods,col=cols,lty=types,lwd=2,horiz=TRUE,bty='n',cex=.88)
 mtext('Log-odds contrast',side=2,outer=TRUE,line=.65,at=.54,cex=1)
}
pdf(file.path(out,'spt_mechanism_36m_2x2.pdf'),width=8.2,height=6.8,useDingbats=FALSE);render();dev.off()
png(file.path(out,'spt_mechanism_36m_2x2.png'),width=2460,height=2040,res=300);render();dev.off()
