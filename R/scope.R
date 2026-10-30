# Exact current main/supplement roster. No catch-all historical method registry.
.table_methods<-c('T-GAM','Source-GAM','Pooled-GAM','Stratified-GAM','SPT','Target-Cox','CoxTL','POTL-Cox',
 'CoxKL','TransCox-BIC','DiscreteKL','Trans-AFT','Random-Survival-Forest','IPCW-XGBoost','Oracle')
paper_scope<-function(design,setting=NULL) {
 if(design %in% c('G2','G3','G5'))return(c(.table_methods,'PR-GAM')) # PR-GAM: DGP-check table only.
 if(design %in% c('level','shape','reverse','structural_D'))return(c('T-GAM','SPT','PR-GAM','Fixed-center','Uncalibrated','SPT-tau0','DST'))
 if(design=='large_n')return(c('T-GAM','SPT-tau0','DST-tau0','SPT','DST'))
 if(design=='clinical')return(c(setdiff(.table_methods,c('Oracle','SPT')),'SPT-GAM','SPT-Cox','SPT-NN','Score-recalibration'))
 stop('Design not displayed in the manuscript: ',design)
}
paper_check_methods<-function(design,methods,setting=NULL) {
 bad<-setdiff(methods,paper_scope(design,setting));if(length(bad))stop('Not displayed for this design: ',paste(bad,collapse=', '))
 stopifnot(!anyDuplicated(methods));invisible(TRUE)
}
