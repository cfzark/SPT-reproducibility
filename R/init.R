# Run entry points from this archive's directory. No parent-project sourcing.
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,TZ="America/Detroit")
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(survival))
source('R/numerics.R');source('R/scope.R');source('R/models.R');source('R/dependencies.R')
.ext_p<-7L;.gcv_trace<-list();.cens_fallbacks<-list()
.e7_helpers<-new.env();.e7_helpers$.make_stratified_folds<-.paper_folds
source('R/simulation.R');source('R/clinical.R')
