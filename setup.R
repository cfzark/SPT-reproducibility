# Explicit dependency setup; no experiment runs and no clinical data downloads.
a<-commandArgs(TRUE);stopifnot(length(a)==1,a%in%c('cran','authors'))
if(a=='cran') {
 install.packages(c('mgcv','survival','jsonlite','digest','timeROC','nnet','ranger','xgboost','Rcpp','RcppArmadillo','RcppEigen','glmnet','mvtnorm','rms','Hmisc','pec','data.table','remotes'),repos='https://cloud.r-project.org')
} else {
 manifest<-jsonlite::fromJSON('config/author_sources.json')
 for(i in seq_len(nrow(manifest))) {
  z<-manifest[i,];path<-file.path('external',z$destination);dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  if(!file.exists(path))download.file(paste0('https://raw.githubusercontent.com/',z$repo,'/',z$commit,'/',z$path),path,mode='wb',method='libcurl')
  stopifnot(digest::digest(file=path,algo='sha256',serialize=FALSE)==z$sha256)
 }
 remotes::install_github('UM-KevinHe/survkl@6d9d776f7f551e192654e4cab3f1a85cc69b3d51',upgrade='never',build_vignettes=FALSE)
 remotes::install_github('UM-KevinHe/DiscreteKL@d742aa71795b6cf0f575a79f4309cd0bbb72d146',upgrade='never',build_vignettes=FALSE)
}
