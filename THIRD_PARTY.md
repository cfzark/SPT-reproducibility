# Comparison implementations

Only methods displayed in the current paper are callable. Author implementations are obtained at pinned commits by `Rscript setup.R authors`; no clinical data, example datasets, or author source archives are bundled.

- CoxTL: YUYING-LU/CoxTL, commit cb2b95b80eda0a2fcc1ba3349a3e03158d42e11a. Density-ratio weighting and native five-fold tuning. The adapter adds tuning metadata to the return object only.
- CoxKL: UM-KevinHe/survkl, commit 6d9d776f7f551e192654e4cab3f1a85cc69b3d51.
- POTL-Cox: yugu-stat/POTL, commit ab9da45d27e5ca37d3445a782ef3d77a67e34d88; application EM.cpp and its three headers. Native profile-likelihood selection.
- DiscreteKL: UM-KevinHe/DiscreteKL, commit d742aa71795b6cf0f575a79f4309cd0bbb72d146.
- Trans-AFT: YuZhengyang-CUEB/Trans-AFT, commit 75b56ca057ee446cbafc6838b7417a240945bd53. Native penalty selection; target-training residual KM gives absolute risks.
- TransCox: local R backend corresponding to Ziyi Li's TransCox 0.1.1, with the historical BIC adapter. R/transcox.R contains only the required R-backend dependency closure. It is not the authors' Python implementation. Original DESCRIPTION and LICENSE are retained verbatim (the former says GPL-2, the supplied license text is GPL-3; this upstream discrepancy remains visible).

Author files keep their upstream licenses. Their download paths and SHA-256 values are in config/author_sources.json. Source retrieval was not repeated during code organization. R/nn.R and the SPT estimator are the project's implementations. This is a local research reproduction archive, not a CRAN package or an online release.
