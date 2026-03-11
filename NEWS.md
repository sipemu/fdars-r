# fdars 0.5.0

## New Features

### Elastic Alignment
- Added `srsf.transform()`, `srsf.inverse()` for Square-Root Slope Function transforms
- Added `elastic.align()` for elastic curve alignment via dynamic programming
- Added `elastic.align.constrained()` for alignment with landmark constraints
- Added `elastic.distance()` for Fisher-Rao distance (self and cross)
- Added `elastic.decomposition()` for amplitude/phase decomposition
- Added `karcher.mean()` for Karcher (Fréchet) mean in the elastic metric
- Added `amplitude.distance()`, `phase.distance()` for component distance matrices
- Added `alignment.quality()` and `alignment.pairwise.consistency()` diagnostics
- Added `periodic.rotate()` with peak, xcorr, landmark, and iterative methods

### Elastic FPCA
- Added `vert.fpca()` — amplitude (vertical) PCA on aligned curves
- Added `horiz.fpca()` — phase (horizontal) PCA on warping functions
- Added `joint.fpca()` — joint amplitude + phase PCA with balance parameter
- S3 print/plot methods for all three classes

### Elastic Regression
- Added `elastic.regression()` — scalar-on-function regression with simultaneous alignment
- Added `elastic.logistic()` — elastic logistic regression for binary classification
- Added `elastic.pcr()` — principal component regression with elastic FPCA
- Added `elastic.attribution()` — amplitude vs phase contribution decomposition
- S3 print/plot methods for all three classes

### Elastic Changepoint Detection
- Added `elastic.changepoint()` — detect structural breaks in functional time series
- Three change types: "amplitude" (shape), "phase" (timing), "fpca" (combined)
- CUSUM test statistics with permutation-based p-values
- S3 print/plot methods with statistic and data visualization

### Model Explainability & Diagnostics
- 22+ explainability methods for `fregre.lm` and `functional.logistic` models:
  - **Global**: `fregre.beta.decomp()`, `fregre.pdp()`, `fregre.ale()`, `fregre.sobol()`, `fregre.importance()`, `fregre.friedman()`
  - **Local**: `fregre.lime()`, `fregre.shap()`, `fregre.counterfactual()`, `fregre.anchor()`, `fregre.saliency()`
  - **Diagnostics**: `fregre.influence()`, `fregre.vif()`, `fregre.dfbetas()`, `fregre.loo()`
  - **Domain**: `fregre.significant.regions()`, `fregre.domain()`, `fregre.pointwise()`
  - **Calibration** (logistic): `fregre.calibration()`, `fregre.ece()`
  - **Uncertainty**: `fregre.prediction.interval()`, `fregre.conformal()`
  - **Meta**: `fregre.stability()`, `fregre.depth()`, `fregre.prototype()`, `fregre.conditional.importance()`

### Penalized Basis Smoothing
- Added `smooth.basis.fd()` — penalized B-spline/Fourier basis smoothing
- Added `smooth.basis.gcv()` — automatic lambda selection via GCV grid search
- S3 print/plot methods for `smooth.basis` objects

### Cross-Validation Framework
- Added `cv.fdata()` — unified k-fold cross-validation for any fitting function
- Supports regression (RMSE, MAE, R²) and classification (accuracy, confusion matrix)
- Repeated k-fold CV with aggregation (mean/SD for regression, majority voting for classification)
- Stratified folding for balanced class splits

### Scalar-on-Function Regression
- Added `fregre.lm()`, `fregre.lm.cv()` — FPC-based functional linear model with Rust backend
- Added `functional.logistic()` — functional logistic regression via IRLS
- Added `fregre.np.mixed()` — mixed scalar + functional nonparametric regression

### Function-on-Scalar Regression
- Added `fosr()` — penalized function-on-scalar regression
- Added `fosr.fpc()` — FPC-based function-on-scalar regression
- Added `fanova()` — functional ANOVA with permutation F-test

### Functional Classification
- Added `fclassif()` — supervised classification (LDA, QDA, kNN, kernel, DD-plot, SVM)
- Added `fclassif.cv()` — cross-validated classification error

### Landmark Registration
- Added `detect.landmarks()` — detect peaks, valleys, zero-crossings, inflection points
- Added `landmark.register()` — register curves by aligning detected landmarks

### TSRVF
- Added `tsrvf.transform()` — tangent-space representation via Karcher mean
- Added `tsrvf.from.alignment()` — build TSRVF from pre-computed alignment
- Added `tsrvf.inverse()` — reconstruct curves from tangent vectors

### Tolerance Bands
- Added `tolerance.band()` with methods: fpca, conformal, scb (Degras), exponential, elastic

### Streaming Depth
- Added `streaming.depth()` — O(T log N) depth computation against reference data

### Soft-DTW Distance
- Added `metric.softDTW()` — differentiable Dynamic Time Warping
- Added `softdtw.barycenter()` — Soft-DTW barycenter computation

### Functional Mixed Models
- Added `fmm()` — functional mixed models with random effects
- Added `fmm.predict()` — prediction from fitted FMM
- Added `fmm.test.fixed()` — permutation test for fixed effects

### GMM Clustering
- Added `cluster.gmm()` — Gaussian Mixture Model clustering with automatic K selection

### Andrews Transformation
- Added `andrews_transform()` — Fourier expansion of multivariate data to functional curves
- Added `andrews_loadings()` — project FPCA eigenfunctions back to original variable space

### Seasonal Analysis
- Period detection: `detect.period()`, `estimate.period()`, `detect.periods()`
- Ensemble methods: `autoperiod()`, `cfd.autoperiod()`, `sazed()`
- Spectral: `lomb.scargle()`, `matrix.profile()`
- Decomposition: `stl.fd()`, `ssa.fd()`, `detrend()`, `decompose()`
- Diagnostics: `seasonal.strength()`, `seasonal.strength.curve()`, `classify.seasonality()`
- Change detection: `detect.seasonality.changes()`, `detect.seasonality.changes.auto()`
- `detect_amplitude_modulation()`, `instantaneous.period()`, `analyze.peak.timing()`

### Simulation Toolbox
- Added `simFunData()`, `simMultiFunData()` — simulate via Karhunen-Loève
- Added `eFun()`, `eVal()` — eigenfunction/eigenvalue generators
- Added `addError()` — add measurement noise
- Added `irregFdata()`, `sparsify()` — irregular/sparse data generation

### Equivalence Testing
- Added `fequiv.test()` — functional equivalence test

### Irregular Functional Data
- Added `irregFdata()` — irregular functional data objects
- Added `as.fdata()` — convert irregular to regular grid

## Improvements
- Upgraded Rust backend (fdars-core) to v0.8.0
- All plots use ggplot2; respects `theme_set()`
- 46 pkgdown articles with SVG overview diagrams
- Unified dispatch: `depth()` (10 methods), `metric()` (14 methods)
- Pre-built binary packages for macOS and Windows (no Rust required)

## Bug Fixes
- Fixed changepoint p-values (normalization mismatch in Brownian bridge simulation)
- Fixed `functional.logistic()` missing `.fpca_scores`, `beta.se`, `std.errors` fields
- Fixed `fdata2basis_2d` example error (extendr wrapper name collision)
- Fixed `mean(fd)` returning matrix instead of fdata
- Fixed TSRVF inverse initial values

## Internal
- Upgraded Rust backend (fdars-core) to v0.8.0
- 200+ new Rust bridge functions
- Vendored Rust crates cleaned for CRAN compliance
- Portable Makefiles (removed GNU extensions)
- GitHub Actions: pinned R version, added Rust toolchain to pkgdown workflow

# fdars 0.4.0

## New Features

### Elastic Alignment
- Added `srsf.transform()` and `srsf.inverse()` for Square-Root Slope Function transforms
- Added `elastic.align()` for elastic curve alignment with S3 print/plot methods
- Added `elastic.distance()` for Fisher-Rao distance computation (self and cross)
- Added `karcher.mean()` for Karcher (Fréchet) mean in the elastic metric with S3 print/plot methods
- Added `metric.elastic()` and "elastic" option to `metric()` dispatcher

### Tolerance Bands
- Added `tolerance.band()` with 5 methods: "fpca", "conformal", "scb", "exponential", "elastic"
- FPCA bootstrap bands (pointwise or simultaneous)
- Distribution-free conformal prediction bands
- Simultaneous confidence bands for the mean (Degras method)
- Exponential family tolerance bands (Gaussian, Binomial, Poisson)
- Elastic tolerance bands (alignment + FPCA)
- S3 print/plot methods for `tolerance.band` objects

### Streaming Depth
- Added `streaming.depth()` for efficient O(T log N) depth computation
- Supports FM, MBD, and BD methods with pre-sorted reference
- Added `depth.streaming()` alias and "streaming" option to `depth()` dispatcher

## Internal
- Upgraded Rust backend (fdars-core) to v0.6.0
- 18 new Rust bridge functions for alignment, tolerance, and streaming depth

# fdars 0.3.2

## Bug Fixes

- Fixed compiled code WARNING: wrapped `abort`/`exit`/`_exit` symbols using linker `--wrap` flag to convert process termination into R errors
- Reduced test CPU time by limiting Rust thread pool (`RAYON_NUM_THREADS=2`)

# fdars 0.3.1

## Bug Fixes

- Fixed Windows installation failure (missing cargo checksum files)
- Wrapped slow bootstrap examples in `\donttest{}`

# fdars 0.3.0

## Internal

- Upgraded Rust backend (fdars-core) to v0.4.0
- New `FdMatrix` type for safer matrix handling (internal)
- New streaming depth module in core (internal)
- Reduced package size by removing non-essential vendored files
- No user-facing API changes

# fdars 0.2.0

## Test Coverage & Quality

- Improved Rust core test coverage to 84%+
- Improved R package test coverage to 80%+
- Added pre-commit hooks for cargo fmt and clippy

## New Features

### Optimal Cluster Selection
- Added `optim.kmeans.fd()` function to automatically determine the optimal number of clusters for functional k-means
- Three selection criteria available:
  - **Silhouette score**: Measures cluster cohesion vs separation (-1 to 1, higher is better)
  - **Calinski-Harabasz index**: Ratio of between/within cluster variance (higher is better)
  - **Elbow method**: Visual inspection of within-cluster sum of squares
- Added `print()` and `plot()` methods for `optim.kmeans.fd` objects
- Silhouette and Calinski-Harabasz computations implemented in Rust for performance

### k-NN Bandwidth Selection for Nonparametric Regression
- Added k-nearest neighbors support to `fregre.np()` via the `type.S` parameter:
  - `"kNN.gCV"`: Global cross-validation (single k for all observations)
  - `"kNN.lCV"`: Local cross-validation (adaptive k per observation)
- Extended `predict.fregre.np()` to handle k-NN models

### Flexible Metrics in Clustering
- `kmeans.fd()` now accepts both string metrics (`"L2"`, `"L1"`, `"Linf"`) and metric/semimetric functions
- String metrics use fast Rust-only path; function metrics provide flexibility for custom distances

## Improvements

### ggplot2 Visualizations
- All plot methods now use ggplot2 instead of base R graphics:
  - `plot.fdata()`: Functional data curves with minimal theme
  - `plot.kmeans.fd()`: Cluster-colored curves with dashed cluster centers
  - `plot.optim.kmeans.fd()`: Criterion scores with optimal k highlighted
  - `plot.outliers.fdata()`: Outlier/normal curves with color legend

## Documentation

### Vignettes
Added 6 comprehensive vignettes:
- Introduction to fdars
- Functional Depth Functions
- Distance Metrics and Semimetrics
- Functional Regression
- Functional Clustering
- Outlier Detection

### API Documentation
- Complete roxygen2 documentation for all exported functions

## Bug Fixes
- Fixed namespace issues with stats and utils imports
- Fixed ggplot2 `.data` pronoun import for R CMD check compliance

# fdars 0.1.0

- Initial release
- Core functional data class (`fdata`) with 1D and 2D support
- 7 depth functions: FM, mode, RP, RT, FSD, KFSD, RPD
- Distance metrics: Lp, Hausdorff, DTW, KL
- Semimetrics: PCA, derivative, basis, Fourier, hshift
- Functional regression: PC, basis, nonparametric
- K-means clustering with k-means++ initialization
- Outlier detection: depth-based and LRT methods
- Statistical tests: flm.test, fmean.test
- Bootstrap inference and confidence intervals
- High-performance Rust backend with parallel processing
