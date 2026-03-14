# Changelog

All notable changes to fdars are documented in this file.

## [0.6.0] - 2026-03-14

### Added

- **Statistical Process Monitoring (SPM)**:
  - `spm.phase1()` — Phase I control chart estimation (FPCA + T²/SPE limits)
  - `spm.monitor()` — Phase II online monitoring with alarm detection
  - `spm.ewma()` — EWMA smoothing for detecting gradual drifts
  - `spm.contributions()` — T²/SPE contribution diagnostics for root cause analysis
  - `mfpca()` — Multivariate FPCA for joint analysis of multiple functional variables
  - `frcc.phase1()`, `frcc.monitor()` — Functional Regression Control Charts
  - S3 print methods for `spm.chart`, `spm.monitor`, `mfpca`, `frcc.chart`

- **Scalar-on-shape regression** (phase-invariant functional regression):
  - `scalar.on.shape()` — regression using elastic alignment and index functions
  - `predict.scalar.on.shape()` — prediction for new curves
  - Index methods: identity, polynomial, Nadaraya-Watson
  - S3 print method for `scalar.on.shape`

- **Robust regression** (outlier-resistant functional regression):
  - `fregre.l1()` — L1 (least absolute deviations) regression
  - `fregre.huber()` — Huber M-estimation with tuning parameter k
  - `predict.fregre.robust()` — prediction from robust models
  - S3 print method for `fregre.robust`

- **Tolerance band extensions**:
  - `tolerance.band(method = "phase")` — tolerance bands for warping functions
  - `tolerance.band(method = "elastic.config")` — joint amplitude + phase elastic bands

- **Calibration indices for conformal prediction**:
  - `conformal.generic.regression()` gains `calibration.indices` parameter
  - `conformal.generic.classification()` gains `calibration.indices` parameter
  - Fixes data leakage when using pre-trained models (GH issue #20)

- **New articles with SVG diagrams** (3 new):
  - Statistical Process Monitoring — Phase I/II control charts, EWMA, FRCC
  - Scalar-on-Shape Regression — elastic alignment + index function regression
  - Robust Regression — L1 and Huber M-estimation under contamination
  - Updated tolerance bands article with phase and elastic.config sections

### Changed
- Upgraded Rust backend (fdars-core) from v0.8.4 to v0.8.5
- `depth.RPD()` now uses Rust backend (was pure R)
- `semimetric.pca()`, `semimetric.deriv()`, `semimetric.basis()` now use Rust backends

### Fixed
- Conformal generic data leakage (#20) — `calibration_indices` parameter enables proper train/calibration split with pre-trained models

## [0.5.0] - 2026-03-11

### Added

- **Elastic alignment** (Fisher-Rao framework):
  - `srsf.transform()`, `srsf.inverse()` — Square-Root Slope Function
  - `elastic.align()` — align curves via dynamic programming
  - `elastic.align.constrained()` — alignment with landmark constraints
  - `elastic.distance()` — Fisher-Rao distance between curves
  - `elastic.decomposition()` — amplitude/phase decomposition
  - `karcher.mean()` — Fréchet mean in the elastic metric
  - `amplitude.distance()`, `phase.distance()` — component distance matrices
  - `alignment.quality()` — diagnostic metrics for alignment results
  - `alignment.pairwise.consistency()` — transitivity checks
  - `periodic.rotate()` — rotate periodic curves with methods: peak, xcorr, landmark, iterative

- **Elastic FPCA** (amplitude/phase PCA):
  - `vert.fpca()` — amplitude (vertical) PCA on aligned curves
  - `horiz.fpca()` — phase (horizontal) PCA on warping functions
  - `joint.fpca()` — joint amplitude + phase PCA with balance parameter
  - S3 print/plot methods for all three classes

- **Elastic regression** (alignment-aware models):
  - `elastic.regression()` — scalar-on-function regression with simultaneous alignment
  - `elastic.logistic()` — elastic logistic regression for binary classification
  - `elastic.pcr()` — principal component regression with elastic FPCA
  - `elastic.attribution()` — amplitude vs phase contribution decomposition
  - S3 print/plot methods for all three classes

- **Elastic changepoint detection**:
  - `elastic.changepoint()` — detect structural breaks in functional time series
  - Three change types: "amplitude" (shape), "phase" (timing), "fpca" (combined)
  - CUSUM test statistics with permutation-based p-values
  - S3 print/plot methods with statistic and data visualization

- **Model explainability & diagnostics** (22+ methods for `fregre.lm` and `functional.logistic`):
  - Global: `fregre.beta.decomp()`, `fregre.pdp()`, `fregre.ale()`, `fregre.sobol()`, `fregre.importance()`, `fregre.friedman()`
  - Local: `fregre.lime()`, `fregre.shap()`, `fregre.counterfactual()`, `fregre.anchor()`, `fregre.saliency()`
  - Diagnostics: `fregre.influence()`, `fregre.vif()`, `fregre.dfbetas()`, `fregre.loo()`
  - Domain: `fregre.significant.regions()`, `fregre.domain()`, `fregre.pointwise()`
  - Calibration (logistic): `fregre.calibration()`, `fregre.ece()`
  - Uncertainty: `fregre.prediction.interval()`, `fregre.conformal()`
  - Meta: `fregre.stability()`, `fregre.depth()`, `fregre.prototype()`, `fregre.conditional.importance()`

- **Penalized basis smoothing**:
  - `smooth.basis.fd()` — penalized B-spline/Fourier basis smoothing
  - `smooth.basis.gcv()` — automatic lambda selection via GCV grid search
  - S3 print/plot methods for `smooth.basis` objects

- **Cross-validation framework**:
  - `cv.fdata()` — unified k-fold cross-validation for any fitting function
  - Supports regression (RMSE, MAE, R²) and classification (accuracy, confusion matrix)
  - Repeated k-fold CV with aggregation
  - Stratified folding for balanced class splits

- **Landmark registration**:
  - `detect.landmarks()` — detect peaks, valleys, zero-crossings, inflection points
  - `landmark.register()` — register curves by aligning detected landmarks

- **TSRVF** (Transported Square-Root Velocity Function):
  - `tsrvf.transform()` — compute tangent-space representation via Karcher mean
  - `tsrvf.from.alignment()` — build TSRVF from pre-computed alignment
  - `tsrvf.inverse()` — reconstruct curves from tangent vectors

- **Tolerance bands**:
  - `tolerance.band()` — functional tolerance/prediction bands with methods: fpca, conformal, scb (Degras), exponential, elastic

- **Streaming depth**:
  - `streaming.depth()` — fast depth computation against pre-sorted reference data

- **Soft-DTW distance**:
  - `metric.softDTW()` — differentiable Dynamic Time Warping
  - `softdtw.barycenter()` — Soft-DTW barycenter computation

- **Function-on-scalar regression**:
  - `fosr()` — penalized function-on-scalar regression
  - `fosr.fpc()` — FPC-based function-on-scalar regression

- **Functional ANOVA**:
  - `fanova()` — one-way functional ANOVA with permutation testing

- **Functional classification**:
  - `fclassif()` — supervised classification (LDA, QDA, kNN, kernel, DD-plot, SVM)
  - `fclassif.cv()` — cross-validated classification error

- **Functional mixed models**:
  - `fmm()` — functional mixed models with subject-level random effects
  - `fmm.predict()` — prediction from fitted FMM
  - `fmm.test.fixed()` — permutation test for fixed effects

- **GMM clustering**:
  - `cluster.gmm()` — Gaussian Mixture Model clustering with automatic K selection
  - `cluster.init()` — k-means++ initialization

- **Scalar-on-function regression**:
  - `fregre.lm()`, `fregre.lm.cv()` — FPC-based functional linear regression with Rust backend
  - `functional.logistic()` — functional logistic regression via IRLS
  - `fregre.np.mixed()` — mixed scalar + functional nonparametric regression

- **Andrews transformation**:
  - `andrews_transform()` — Fourier expansion of multivariate data to functional curves
  - `andrews_loadings()` — project FPCA eigenfunctions back to original variable space

- **Seasonal analysis** (comprehensive module):
  - Period detection: `detect.period()`, `estimate.period()`, `detect.periods()`
  - Ensemble methods: `autoperiod()`, `cfd.autoperiod()`, `sazed()`
  - Spectral: `lomb.scargle()`, `matrix.profile()`
  - Decomposition: `stl.fd()`, `ssa.fd()`, `detrend()`, `decompose()`
  - Diagnostics: `seasonal.strength()`, `seasonal.strength.curve()`, `classify.seasonality()`
  - Change detection: `detect.seasonality.changes()`, `detect.seasonality.changes.auto()`
  - `detect_amplitude_modulation()`, `instantaneous.period()`, `analyze.peak.timing()`
  - S3 plot methods: `plot.peak_detection()`, `plot.peak_timing()`

- **Equivalence testing**:
  - `fequiv.test()` — functional equivalence test

- **Simulation toolbox**:
  - `simFunData()`, `simMultiFunData()` — simulate functional data via Karhunen-Loève
  - `eFun()`, `eVal()` — eigenfunction/eigenvalue generators
  - `addError()` — add measurement noise
  - `irregFdata()`, `sparsify()` — irregular/sparse functional data

- **Irregular functional data**:
  - `irregFdata()` — irregular functional data objects
  - `as.fdata()` — convert irregular to regular grid
  - `mean.irregFdata()`, `summary.irregFdata()`, `plot.irregFdata()`

- **Unified dispatch interfaces**:
  - `depth()` — dispatches to FM, mode, RP, RT, BD, MBD, MEI, FSD, KFSD, RPD
  - `metric()` — dispatches to lp, hausdorff, dtw, softdtw, kl, elastic, amplitude, phase, pca, deriv, basis, fourier, hshift

- **Plot methods** for all new object types:
  - `plot.elastic.align()`, `plot.karcher.mean()`, `plot.alignment.quality()`
  - `plot.landmark.register()`, `plot.tsrvf()`
  - `plot.tolerance.band()`
  - `plot.fclassif()`, `plot.fmm()`, `plot.fosr()`, `plot.fanova()`
  - `plot.cluster.gmm()`
  - `plot.elastic.changepoint()`, `plot.elastic.regression()`, `plot.elastic.logistic()`, `plot.elastic.pcr()`
  - `plot.vert.fpca()`, `plot.horiz.fpca()`, `plot.joint.fpca()`
  - `plot.smooth.basis()`, `plot.cv.fdata()`

- **Vignettes/articles** (46 articles with SVG overview diagrams):
  - Elastic alignment, landmark registration, TSRVF, alignment comparison
  - Elastic FPCA, elastic regression, elastic changepoint detection
  - Tolerance bands, streaming depth, distance metrics
  - Functional regression, classification, mixed models, GMM clustering
  - Model explainability, regression diagnostics, uncertainty quantification
  - Penalized basis smoothing, cross-validation framework
  - Andrews transformation, seasonal analysis
  - Covariance functions, simulation toolbox, custom plotting
  - Example analyses: Berkeley growth study, wine dataset, Sonar TSRVF classification

- **Pre-built binary packages** for macOS (.tgz) and Windows (.zip) attached to
  GitHub Releases — install without Rust toolchain

### Changed
- Upgraded Rust backend (fdars-core) to v0.8.0
- All plots use ggplot2; respects `theme_set()`
- Renamed functions to avoid S3 dispatch conflicts:
  - `fdata2basis.cv()` → `fdata2basis_cv()`
  - `fdata2basis.2d()` → `fdata2basis_2d()`
  - `basis2fdata.2d()` → `basis2fdata_2d()`
  - `norm.fdata()` → `norm()`
  - `cov.*()` → `kernel.*()` (covariance kernel functions)
- Internal Rust wrapper functions no longer generate man pages (CRAN compliance)
- Portable Makefiles: removed GNU extensions (`export`, `:=`)
- Vendored Rust crates cleaned of hidden files and long paths
- GitHub Actions: pinned R 4.5.2, added Rust toolchain to pkgdown workflow

### Fixed
- Changepoint p-values always 1.0 (normalization mismatch in Brownian bridge simulation; replaced with permutation testing)
- `functional.logistic()` missing `.fpca_scores`, `beta.se`, `std.errors` from Rust
- `fregre.conformal()` wrapper now works correctly
- `fdata2basis_2d` example error caused by extendr wrapper name collision
- `mean(fd)` now returns fdata object (was returning matrix)
- TSRVF inverse initial values
- Empty bar charts in seasonal plots (use `coord_cartesian` instead of `ylim`)
- Regression article cross-references and See Also sections

## [0.4.0] - 2026-03-04

### Added
- `id` and `metadata` slots in fdata objects for storing curve identifiers and associated data
- Outlier plot labeling: `plot(outliergram, label = "id")` or `label = "column_name"`
- `magnitudeshape()` labeling support (renamed from MS.plot)
- `fregre.np.multi()` for regression with multiple functional predictors
- Enhanced `plot.fdata()` with group coloring, mean curves, and confidence intervals
- `group.distance()` for measuring distances between groups of curves
- `group.test()` permutation test for group differences
- `outliergram()` visualization (MEI vs MBD plot)
- `plot.fdata2pc()` for FPCA visualization (components, variance, scores)

### Changed
- Auto-reduce alpha when `show.mean = TRUE` in `plot.fdata()`
- Renamed `fdata.deriv()` to `deriv()` for consistency

### Fixed
- `plot.group.distance()` error handling
- Release workflow now generates documentation before building
- Missing `%||%` operator definition

## [0.3.0] - 2026-02-19

### Added
- Covariance kernel functions: `kernel.exponential()`, `kernel.matern()`, `kernel.gaussian()`, etc.
- `make.gaussian.process()` for simulating Gaussian process realizations
- 2D functional data support for most functions
- Unified API: `depth()`, `median()`, `trimmed()`, `trimvar()` with method parameter
- Band depth (`depth.BD`, `depth.MBD`, `depth.MEI`)
- Functional boxplot (`boxplot.fdata`)
- MS-plot for outlier detection
- Fuzzy c-means clustering (`cluster.fcm`)
- Geometric median (`gmed`)
- Curve registration (`register.fd`)
- Local averages feature extraction (`localavg.fdata`)
- Optimal k selection for k-means (`cluster.optim`)
- k-NN bandwidth selection for nonparametric regression
- **Basis representation module** with Rust backend:
  - `fdata2basis()`, `basis2fdata()` — convert between data and basis coefficients
  - `basis.gcv()`, `basis.aic()`, `basis.bic()` — goodness-of-fit metrics
  - `fdata2basis_cv()` — cross-validation for optimal nbasis selection
  - `pspline()` — P-spline smoothing with automatic lambda selection
  - `fdata2basis_2d()`, `basis2fdata_2d()` — 2D tensor product basis support
  - `pspline.2d()` — 2D P-spline with anisotropic penalties
  - `select.basis.auto()` — automatic basis selection per curve
  - `normalize()` — scale curves to unit Lp norm

### Changed
- Cleaned up API: removed backward compatibility shims
- Renamed functions for consistency (e.g., `fdata.mean` → `mean.fdata`)
- All plots now use ggplot2 instead of base R graphics

## [0.1.0] - 2026-02-13

### Initial Release
- **Core functional data structure**: `fdata()` for 1D and 2D functional data
- **Depth functions**: FM, mode, RP, RT, FSD, KFSD, RPD (all with Rust backend)
- **Statistics**: mean, variance, standard deviation, covariance
- **Distance metrics**: Lp, Hausdorff, DTW, KL divergence
- **Semimetrics**: basis projection, Fourier, horizontal shift, PCA, derivative
- **Regression**: `fregre.pc()`, `fregre.basis()`, `fregre.np()` with CV variants
- **Outlier detection**: `outliers.depth.pond()`, `outliers.depth.trim()`, `outliers.lrt()`
- **Smoothing**: Nadaraya-Watson, local linear, local polynomial, k-NN
- **Clustering**: functional k-means
- **Hypothesis testing**: `flm.test()`, `fmean.test.fdata()`
- **Utilities**: Simpson integration, inner product, derivatives
