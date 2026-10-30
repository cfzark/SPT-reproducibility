# SPT paper: simulation and clinical code

This directory contains only analyses in **main.tex and supplement.tex as of 2026-10-29**. It is independent of the historical project `src/` tree. No patient data, fitted clinical models, unused transfer variants, old cohorts, or historical experiment archives are included. The general-purpose R package is separate at `packages/SPT` in the project.

## Included scope

| Entry | Settings | Methods |
|---|---|---|
| G2 main simulation | A–D, log-logistic/Weibull, nT=500 | 15 methods in the full supplement tables; PR-GAM additionally supplies the DGP-check table |
| G3 sample size | B,D; nT=200/1000 | Same displayed roster |
| G5 irrelevant predictors | B,D; p=20, nT=500 | Same displayed roster |
| Source mismatch | level, shape, reverse | T-GAM and the six displayed component comparisons |
| Structural mismatch | GCV→REML D only | Same six component comparisons, including DST |
| Large-n figure | fixed-additive log-logistic D; nT=500/1000/2000/4000/10000 | T-GAM, SPT/DST at tau=0 and one-SE-half |
| Clinical | UKB B → SPRINT standard; death; seven predictors; 24/30/36 months | Exactly 17 methods in the current supplementary clinical table |

## Quick use

```sh
Rscript setup.R cran       # install required CRAN packages, if needed
Rscript setup.R authors    # optional pinned comparison implementations
Rscript tests/smoke.R

# One complete main-design replication (cell IDs are in config/simulation_cells.csv):
Rscript run_simulation.R G2 1 1 /your/output/G2_A_LL_rep001 all

# Same setting, only displayed GAM methods:
Rscript run_simulation.R G2 1 1 /your/output/G2_A_LL_core001 T-GAM,SPT,PR-GAM,Source-GAM,Pooled-GAM,Stratified-GAM,Oracle

# Separate D matched ablation; fixed-additive large-n D:
Rscript run_simulation.R structural_D 4 1 /your/output/D_ablation001 all
SPT_SAVE_MODELS=1 Rscript run_simulation.R large_n 4 1 /your/output/D_n2000_rep001 all 2000
Rscript diagnose_large_n.R /your/output/D_n2000_rep001/simulation_models.rds /your/output/diagnostics001.csv

# Clinical: prepare shared sources once, then run repeats 1,...,100 independently.
Rscript prepare_realdata_sources.R /your/cleaned_data /your/private/source_cache
SPT_SAVE_MODELS=1 Rscript run_realdata.R /your/cleaned_data 1 /your/private/rep_001 /your/private/source_cache all
```

Replace `rep=1` by 1–100 for the paper design. Source mismatch uses cell 1 (log-logistic) or 5 (Weibull) with `level`, `shape` or `reverse`. G3/G5 cell IDs are explicit in the calibration CSV. Runs reject existing metric files rather than overwrite them. Selected-method jobs can be separated into different output directories. Preparing the clinical source cache serially avoids duplicate NN tuning.
