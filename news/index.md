# Changelog

## fdars 0.5.0

### New Features

#### Elastic Alignment

- Added
  [`srsf.transform()`](https://sipemu.github.io/fdars-r/reference/srsf.transform.md),
  [`srsf.inverse()`](https://sipemu.github.io/fdars-r/reference/srsf.inverse.md)
  for Square-Root Slope Function transforms
- Added
  [`elastic.align()`](https://sipemu.github.io/fdars-r/reference/elastic.align.md)
  for elastic curve alignment via dynamic programming
- Added
  [`elastic.align.constrained()`](https://sipemu.github.io/fdars-r/reference/elastic.align.constrained.md)
  for alignment with landmark constraints
- Added
  [`elastic.distance()`](https://sipemu.github.io/fdars-r/reference/elastic.distance.md)
  for Fisher-Rao distance (self and cross)
- Added
  [`elastic.decomposition()`](https://sipemu.github.io/fdars-r/reference/elastic.decomposition.md)
  for amplitude/phase decomposition
- Added
  [`karcher.mean()`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md)
  for Karcher (Fréchet) mean in the elastic metric
- Added
  [`amplitude.distance()`](https://sipemu.github.io/fdars-r/reference/amplitude.distance.md),
  [`phase.distance()`](https://sipemu.github.io/fdars-r/reference/phase.distance.md)
  for component distance matrices
- Added
  [`alignment.quality()`](https://sipemu.github.io/fdars-r/reference/alignment.quality.md)
  and
  [`alignment.pairwise.consistency()`](https://sipemu.github.io/fdars-r/reference/alignment.pairwise.consistency.md)
  diagnostics
- Added
  [`periodic.rotate()`](https://sipemu.github.io/fdars-r/reference/periodic.rotate.md)
  with peak, xcorr, landmark, and iterative methods

#### Elastic FPCA

- Added
  [`vert.fpca()`](https://sipemu.github.io/fdars-r/reference/vert.fpca.md)
  — amplitude (vertical) PCA on aligned curves
- Added
  [`horiz.fpca()`](https://sipemu.github.io/fdars-r/reference/horiz.fpca.md)
  — phase (horizontal) PCA on warping functions
- Added
  [`joint.fpca()`](https://sipemu.github.io/fdars-r/reference/joint.fpca.md)
  — joint amplitude + phase PCA with balance parameter
- S3 print/plot methods for all three classes

#### Elastic Regression

- Added
  [`elastic.regression()`](https://sipemu.github.io/fdars-r/reference/elastic.regression.md)
  — scalar-on-function regression with simultaneous alignment
- Added
  [`elastic.logistic()`](https://sipemu.github.io/fdars-r/reference/elastic.logistic.md)
  — elastic logistic regression for binary classification
- Added
  [`elastic.pcr()`](https://sipemu.github.io/fdars-r/reference/elastic.pcr.md)
  — principal component regression with elastic FPCA
- Added
  [`elastic.attribution()`](https://sipemu.github.io/fdars-r/reference/elastic.attribution.md)
  — amplitude vs phase contribution decomposition
- S3 print/plot methods for all three classes

#### Elastic Changepoint Detection

- Added
  [`elastic.changepoint()`](https://sipemu.github.io/fdars-r/reference/elastic.changepoint.md)
  — detect structural breaks in functional time series
- Three change types: “amplitude” (shape), “phase” (timing), “fpca”
  (combined)
- CUSUM test statistics with permutation-based p-values
- S3 print/plot methods with statistic and data visualization

#### Model Explainability & Diagnostics

- 22+ explainability methods for `fregre.lm` and `functional.logistic`
  models:
  - **Global**:
    [`fregre.beta.decomp()`](https://sipemu.github.io/fdars-r/reference/fregre.beta.decomp.md),
    [`fregre.pdp()`](https://sipemu.github.io/fdars-r/reference/fregre.pdp.md),
    [`fregre.ale()`](https://sipemu.github.io/fdars-r/reference/fregre.ale.md),
    [`fregre.sobol()`](https://sipemu.github.io/fdars-r/reference/fregre.sobol.md),
    [`fregre.importance()`](https://sipemu.github.io/fdars-r/reference/fregre.importance.md),
    [`fregre.friedman()`](https://sipemu.github.io/fdars-r/reference/fregre.friedman.md)
  - **Local**:
    [`fregre.lime()`](https://sipemu.github.io/fdars-r/reference/fregre.lime.md),
    [`fregre.shap()`](https://sipemu.github.io/fdars-r/reference/fregre.shap.md),
    [`fregre.counterfactual()`](https://sipemu.github.io/fdars-r/reference/fregre.counterfactual.md),
    [`fregre.anchor()`](https://sipemu.github.io/fdars-r/reference/fregre.anchor.md),
    [`fregre.saliency()`](https://sipemu.github.io/fdars-r/reference/fregre.saliency.md)
  - **Diagnostics**:
    [`fregre.influence()`](https://sipemu.github.io/fdars-r/reference/fregre.influence.md),
    [`fregre.vif()`](https://sipemu.github.io/fdars-r/reference/fregre.vif.md),
    [`fregre.dfbetas()`](https://sipemu.github.io/fdars-r/reference/fregre.dfbetas.md),
    [`fregre.loo()`](https://sipemu.github.io/fdars-r/reference/fregre.loo.md)
  - **Domain**:
    [`fregre.significant.regions()`](https://sipemu.github.io/fdars-r/reference/fregre.significant.regions.md),
    [`fregre.domain()`](https://sipemu.github.io/fdars-r/reference/fregre.domain.md),
    [`fregre.pointwise()`](https://sipemu.github.io/fdars-r/reference/fregre.pointwise.md)
  - **Calibration** (logistic):
    [`fregre.calibration()`](https://sipemu.github.io/fdars-r/reference/fregre.calibration.md),
    [`fregre.ece()`](https://sipemu.github.io/fdars-r/reference/fregre.ece.md)
  - **Uncertainty**:
    [`fregre.prediction.interval()`](https://sipemu.github.io/fdars-r/reference/fregre.prediction.interval.md),
    [`fregre.conformal()`](https://sipemu.github.io/fdars-r/reference/fregre.conformal.md)
  - **Meta**:
    [`fregre.stability()`](https://sipemu.github.io/fdars-r/reference/fregre.stability.md),
    [`fregre.depth()`](https://sipemu.github.io/fdars-r/reference/fregre.depth.md),
    [`fregre.prototype()`](https://sipemu.github.io/fdars-r/reference/fregre.prototype.md),
    [`fregre.conditional.importance()`](https://sipemu.github.io/fdars-r/reference/fregre.conditional.importance.md)

#### Penalized Basis Smoothing

- Added
  [`smooth.basis.fd()`](https://sipemu.github.io/fdars-r/reference/smooth.basis.fd.md)
  — penalized B-spline/Fourier basis smoothing
- Added
  [`smooth.basis.gcv()`](https://sipemu.github.io/fdars-r/reference/smooth.basis.gcv.md)
  — automatic lambda selection via GCV grid search
- S3 print/plot methods for `smooth.basis` objects

#### Cross-Validation Framework

- Added
  [`cv.fdata()`](https://sipemu.github.io/fdars-r/reference/cv.fdata.md)
  — unified k-fold cross-validation for any fitting function
- Supports regression (RMSE, MAE, R²) and classification (accuracy,
  confusion matrix)
- Repeated k-fold CV with aggregation (mean/SD for regression, majority
  voting for classification)
- Stratified folding for balanced class splits

#### Scalar-on-Function Regression

- Added
  [`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md),
  [`fregre.lm.cv()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.cv.md)
  — FPC-based functional linear model with Rust backend
- Added
  [`functional.logistic()`](https://sipemu.github.io/fdars-r/reference/functional.logistic.md)
  — functional logistic regression via IRLS
- Added
  [`fregre.np.mixed()`](https://sipemu.github.io/fdars-r/reference/fregre.np.mixed.md)
  — mixed scalar + functional nonparametric regression

#### Function-on-Scalar Regression

- Added [`fosr()`](https://sipemu.github.io/fdars-r/reference/fosr.md) —
  penalized function-on-scalar regression
- Added
  [`fosr.fpc()`](https://sipemu.github.io/fdars-r/reference/fosr.fpc.md)
  — FPC-based function-on-scalar regression
- Added
  [`fanova()`](https://sipemu.github.io/fdars-r/reference/fanova.md) —
  functional ANOVA with permutation F-test

#### Functional Classification

- Added
  [`fclassif()`](https://sipemu.github.io/fdars-r/reference/fclassif.md)
  — supervised classification (LDA, QDA, kNN, kernel, DD-plot, SVM)
- Added
  [`fclassif.cv()`](https://sipemu.github.io/fdars-r/reference/fclassif.cv.md)
  — cross-validated classification error

#### Landmark Registration

- Added
  [`detect.landmarks()`](https://sipemu.github.io/fdars-r/reference/detect.landmarks.md)
  — detect peaks, valleys, zero-crossings, inflection points
- Added
  [`landmark.register()`](https://sipemu.github.io/fdars-r/reference/landmark.register.md)
  — register curves by aligning detected landmarks

#### TSRVF

- Added
  [`tsrvf.transform()`](https://sipemu.github.io/fdars-r/reference/tsrvf.transform.md)
  — tangent-space representation via Karcher mean
- Added
  [`tsrvf.from.alignment()`](https://sipemu.github.io/fdars-r/reference/tsrvf.from.alignment.md)
  — build TSRVF from pre-computed alignment
- Added
  [`tsrvf.inverse()`](https://sipemu.github.io/fdars-r/reference/tsrvf.inverse.md)
  — reconstruct curves from tangent vectors

#### Tolerance Bands

- Added
  [`tolerance.band()`](https://sipemu.github.io/fdars-r/reference/tolerance.band.md)
  with methods: fpca, conformal, scb (Degras), exponential, elastic

#### Streaming Depth

- Added
  [`streaming.depth()`](https://sipemu.github.io/fdars-r/reference/streaming.depth.md)
  — O(T log N) depth computation against reference data

#### Soft-DTW Distance

- Added
  [`metric.softDTW()`](https://sipemu.github.io/fdars-r/reference/metric.softDTW.md)
  — differentiable Dynamic Time Warping
- Added
  [`softdtw.barycenter()`](https://sipemu.github.io/fdars-r/reference/softdtw.barycenter.md)
  — Soft-DTW barycenter computation

#### Functional Mixed Models

- Added [`fmm()`](https://sipemu.github.io/fdars-r/reference/fmm.md) —
  functional mixed models with random effects
- Added
  [`fmm.predict()`](https://sipemu.github.io/fdars-r/reference/fmm.predict.md)
  — prediction from fitted FMM
- Added
  [`fmm.test.fixed()`](https://sipemu.github.io/fdars-r/reference/fmm.test.fixed.md)
  — permutation test for fixed effects

#### GMM Clustering

- Added
  [`cluster.gmm()`](https://sipemu.github.io/fdars-r/reference/cluster.gmm.md)
  — Gaussian Mixture Model clustering with automatic K selection

#### Andrews Transformation

- Added
  [`andrews_transform()`](https://sipemu.github.io/fdars-r/reference/andrews_transform.md)
  — Fourier expansion of multivariate data to functional curves
- Added
  [`andrews_loadings()`](https://sipemu.github.io/fdars-r/reference/andrews_loadings.md)
  — project FPCA eigenfunctions back to original variable space

#### Seasonal Analysis

- Period detection:
  [`detect.period()`](https://sipemu.github.io/fdars-r/reference/detect.period.md),
  [`estimate.period()`](https://sipemu.github.io/fdars-r/reference/estimate.period.md),
  [`detect.periods()`](https://sipemu.github.io/fdars-r/reference/detect.periods.md)
- Ensemble methods:
  [`autoperiod()`](https://sipemu.github.io/fdars-r/reference/autoperiod.md),
  [`cfd.autoperiod()`](https://sipemu.github.io/fdars-r/reference/cfd.autoperiod.md),
  [`sazed()`](https://sipemu.github.io/fdars-r/reference/sazed.md)
- Spectral:
  [`lomb.scargle()`](https://sipemu.github.io/fdars-r/reference/lomb.scargle.md),
  [`matrix.profile()`](https://sipemu.github.io/fdars-r/reference/matrix.profile.md)
- Decomposition:
  [`stl.fd()`](https://sipemu.github.io/fdars-r/reference/stl.fd.md),
  [`ssa.fd()`](https://sipemu.github.io/fdars-r/reference/ssa.fd.md),
  [`detrend()`](https://sipemu.github.io/fdars-r/reference/detrend.md),
  [`decompose()`](https://sipemu.github.io/fdars-r/reference/decompose.md)
- Diagnostics:
  [`seasonal.strength()`](https://sipemu.github.io/fdars-r/reference/seasonal.strength.md),
  [`seasonal.strength.curve()`](https://sipemu.github.io/fdars-r/reference/seasonal.strength.curve.md),
  [`classify.seasonality()`](https://sipemu.github.io/fdars-r/reference/classify.seasonality.md)
- Change detection:
  [`detect.seasonality.changes()`](https://sipemu.github.io/fdars-r/reference/detect.seasonality.changes.md),
  [`detect.seasonality.changes.auto()`](https://sipemu.github.io/fdars-r/reference/detect.seasonality.changes.auto.md)
- [`detect_amplitude_modulation()`](https://sipemu.github.io/fdars-r/reference/detect_amplitude_modulation.md),
  [`instantaneous.period()`](https://sipemu.github.io/fdars-r/reference/instantaneous.period.md),
  [`analyze.peak.timing()`](https://sipemu.github.io/fdars-r/reference/analyze.peak.timing.md)

#### Simulation Toolbox

- Added
  [`simFunData()`](https://sipemu.github.io/fdars-r/reference/simFunData.md),
  [`simMultiFunData()`](https://sipemu.github.io/fdars-r/reference/simMultiFunData.md)
  — simulate via Karhunen-Loève
- Added [`eFun()`](https://sipemu.github.io/fdars-r/reference/eFun.md),
  [`eVal()`](https://sipemu.github.io/fdars-r/reference/eVal.md) —
  eigenfunction/eigenvalue generators
- Added
  [`addError()`](https://sipemu.github.io/fdars-r/reference/addError.md)
  — add measurement noise
- Added
  [`irregFdata()`](https://sipemu.github.io/fdars-r/reference/irregFdata.md),
  [`sparsify()`](https://sipemu.github.io/fdars-r/reference/sparsify.md)
  — irregular/sparse data generation

#### Equivalence Testing

- Added
  [`fequiv.test()`](https://sipemu.github.io/fdars-r/reference/fequiv.test.md)
  — functional equivalence test

#### Irregular Functional Data

- Added
  [`irregFdata()`](https://sipemu.github.io/fdars-r/reference/irregFdata.md)
  — irregular functional data objects
- Added
  [`as.fdata()`](https://sipemu.github.io/fdars-r/reference/as.fdata.irregFdata.md)
  — convert irregular to regular grid

### Improvements

- Upgraded Rust backend (fdars-core) to v0.8.0
- All plots use ggplot2; respects
  [`theme_set()`](https://ggplot2.tidyverse.org/reference/get_theme.html)
- 46 pkgdown articles with SVG overview diagrams
- Unified dispatch:
  [`depth()`](https://sipemu.github.io/fdars-r/reference/depth.md) (10
  methods),
  [`metric()`](https://sipemu.github.io/fdars-r/reference/metric.md) (14
  methods)
- Pre-built binary packages for macOS and Windows (no Rust required)

### Bug Fixes

- Fixed changepoint p-values (normalization mismatch in Brownian bridge
  simulation)
- Fixed
  [`functional.logistic()`](https://sipemu.github.io/fdars-r/reference/functional.logistic.md)
  missing `.fpca_scores`, `beta.se`, `std.errors` fields
- Fixed `fdata2basis_2d` example error (extendr wrapper name collision)
- Fixed `mean(fd)` returning matrix instead of fdata
- Fixed TSRVF inverse initial values

### Internal

- Upgraded Rust backend (fdars-core) to v0.8.0
- 200+ new Rust bridge functions
- Vendored Rust crates cleaned for CRAN compliance
- Portable Makefiles (removed GNU extensions)
- GitHub Actions: pinned R version, added Rust toolchain to pkgdown
  workflow

## fdars 0.4.0

### New Features

#### Elastic Alignment

- Added
  [`srsf.transform()`](https://sipemu.github.io/fdars-r/reference/srsf.transform.md)
  and
  [`srsf.inverse()`](https://sipemu.github.io/fdars-r/reference/srsf.inverse.md)
  for Square-Root Slope Function transforms
- Added
  [`elastic.align()`](https://sipemu.github.io/fdars-r/reference/elastic.align.md)
  for elastic curve alignment with S3 print/plot methods
- Added
  [`elastic.distance()`](https://sipemu.github.io/fdars-r/reference/elastic.distance.md)
  for Fisher-Rao distance computation (self and cross)
- Added
  [`karcher.mean()`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md)
  for Karcher (Fréchet) mean in the elastic metric with S3 print/plot
  methods
- Added
  [`metric.elastic()`](https://sipemu.github.io/fdars-r/reference/metric.elastic.md)
  and “elastic” option to
  [`metric()`](https://sipemu.github.io/fdars-r/reference/metric.md)
  dispatcher

#### Tolerance Bands

- Added
  [`tolerance.band()`](https://sipemu.github.io/fdars-r/reference/tolerance.band.md)
  with 5 methods: “fpca”, “conformal”, “scb”, “exponential”, “elastic”
- FPCA bootstrap bands (pointwise or simultaneous)
- Distribution-free conformal prediction bands
- Simultaneous confidence bands for the mean (Degras method)
- Exponential family tolerance bands (Gaussian, Binomial, Poisson)
- Elastic tolerance bands (alignment + FPCA)
- S3 print/plot methods for `tolerance.band` objects

#### Streaming Depth

- Added
  [`streaming.depth()`](https://sipemu.github.io/fdars-r/reference/streaming.depth.md)
  for efficient O(T log N) depth computation
- Supports FM, MBD, and BD methods with pre-sorted reference
- Added
  [`depth.streaming()`](https://sipemu.github.io/fdars-r/reference/depth.streaming.md)
  alias and “streaming” option to
  [`depth()`](https://sipemu.github.io/fdars-r/reference/depth.md)
  dispatcher

### Internal

- Upgraded Rust backend (fdars-core) to v0.6.0
- 18 new Rust bridge functions for alignment, tolerance, and streaming
  depth

## fdars 0.3.2

### Bug Fixes

- Fixed compiled code WARNING: wrapped `abort`/`exit`/`_exit` symbols
  using linker `--wrap` flag to convert process termination into R
  errors
- Reduced test CPU time by limiting Rust thread pool
  (`RAYON_NUM_THREADS=2`)

## fdars 0.3.1

### Bug Fixes

- Fixed Windows installation failure (missing cargo checksum files)
- Wrapped slow bootstrap examples in `\donttest{}`

## fdars 0.3.0

### Internal

- Upgraded Rust backend (fdars-core) to v0.4.0
- New `FdMatrix` type for safer matrix handling (internal)
- New streaming depth module in core (internal)
- Reduced package size by removing non-essential vendored files
- No user-facing API changes

## fdars 0.2.0

### Test Coverage & Quality

- Improved Rust core test coverage to 84%+
- Improved R package test coverage to 80%+
- Added pre-commit hooks for cargo fmt and clippy

### New Features

#### Optimal Cluster Selection

- Added `optim.kmeans.fd()` function to automatically determine the
  optimal number of clusters for functional k-means
- Three selection criteria available:
  - **Silhouette score**: Measures cluster cohesion vs separation (-1 to
    1, higher is better)
  - **Calinski-Harabasz index**: Ratio of between/within cluster
    variance (higher is better)
  - **Elbow method**: Visual inspection of within-cluster sum of squares
- Added [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for
  `optim.kmeans.fd` objects
- Silhouette and Calinski-Harabasz computations implemented in Rust for
  performance

#### k-NN Bandwidth Selection for Nonparametric Regression

- Added k-nearest neighbors support to
  [`fregre.np()`](https://sipemu.github.io/fdars-r/reference/fregre.np.md)
  via the `type.S` parameter:
  - `"kNN.gCV"`: Global cross-validation (single k for all observations)
  - `"kNN.lCV"`: Local cross-validation (adaptive k per observation)
- Extended
  [`predict.fregre.np()`](https://sipemu.github.io/fdars-r/reference/predict.fregre.np.md)
  to handle k-NN models

#### Flexible Metrics in Clustering

- `kmeans.fd()` now accepts both string metrics (`"L2"`, `"L1"`,
  `"Linf"`) and metric/semimetric functions
- String metrics use fast Rust-only path; function metrics provide
  flexibility for custom distances

### Improvements

#### ggplot2 Visualizations

- All plot methods now use ggplot2 instead of base R graphics:
  - [`plot.fdata()`](https://sipemu.github.io/fdars-r/reference/plot.fdata.md):
    Functional data curves with minimal theme
  - `plot.kmeans.fd()`: Cluster-colored curves with dashed cluster
    centers
  - `plot.optim.kmeans.fd()`: Criterion scores with optimal k
    highlighted
  - [`plot.outliers.fdata()`](https://sipemu.github.io/fdars-r/reference/plot.outliers.fdata.md):
    Outlier/normal curves with color legend

### Documentation

#### Vignettes

Added 6 comprehensive vignettes: - Introduction to fdars - Functional
Depth Functions - Distance Metrics and Semimetrics - Functional
Regression - Functional Clustering - Outlier Detection

#### API Documentation

- Complete roxygen2 documentation for all exported functions

### Bug Fixes

- Fixed namespace issues with stats and utils imports
- Fixed ggplot2 `.data` pronoun import for R CMD check compliance

## fdars 0.1.0

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
