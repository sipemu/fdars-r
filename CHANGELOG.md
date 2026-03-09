# Changelog

All notable changes to fdars are documented in this file.

## \[0.5.0\] - 2025-03-09

### Added

- **Elastic alignment** (Fisher-Rao framework):
  - [`srsf.transform()`](https://sipemu.github.io/fdars-r/reference/srsf.transform.md),
    [`srsf.inverse()`](https://sipemu.github.io/fdars-r/reference/srsf.inverse.md)
    — Square-Root Slope Function
  - [`elastic.align()`](https://sipemu.github.io/fdars-r/reference/elastic.align.md)
    — align curves via dynamic programming
  - [`elastic.align.constrained()`](https://sipemu.github.io/fdars-r/reference/elastic.align.constrained.md)
    — alignment with landmark constraints
  - [`elastic.distance()`](https://sipemu.github.io/fdars-r/reference/elastic.distance.md)
    — Fisher-Rao distance between curves
  - [`elastic.decomposition()`](https://sipemu.github.io/fdars-r/reference/elastic.decomposition.md)
    — amplitude/phase decomposition
  - [`karcher.mean()`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md)
    — Fréchet mean in the elastic metric
  - [`amplitude.distance()`](https://sipemu.github.io/fdars-r/reference/amplitude.distance.md),
    [`phase.distance()`](https://sipemu.github.io/fdars-r/reference/phase.distance.md)
    — component distance matrices
  - [`alignment.quality()`](https://sipemu.github.io/fdars-r/reference/alignment.quality.md)
    — diagnostic metrics for alignment results
  - [`alignment.pairwise.consistency()`](https://sipemu.github.io/fdars-r/reference/alignment.pairwise.consistency.md)
    — transitivity checks
  - [`periodic.rotate()`](https://sipemu.github.io/fdars-r/reference/periodic.rotate.md)
    — rotate periodic curves with methods: peak, xcorr, landmark,
    iterative
- **Landmark registration**:
  - [`detect.landmarks()`](https://sipemu.github.io/fdars-r/reference/detect.landmarks.md)
    — detect peaks, valleys, zero-crossings, inflection points
  - [`landmark.register()`](https://sipemu.github.io/fdars-r/reference/landmark.register.md)
    — register curves by aligning detected landmarks
- **TSRVF** (Transported Square-Root Velocity Function):
  - [`tsrvf.transform()`](https://sipemu.github.io/fdars-r/reference/tsrvf.transform.md)
    — compute tangent-space representation via Karcher mean
  - [`tsrvf.from.alignment()`](https://sipemu.github.io/fdars-r/reference/tsrvf.from.alignment.md)
    — build TSRVF from pre-computed alignment
  - [`tsrvf.inverse()`](https://sipemu.github.io/fdars-r/reference/tsrvf.inverse.md)
    — reconstruct curves from tangent vectors
- **Tolerance bands**:
  - [`tolerance.band()`](https://sipemu.github.io/fdars-r/reference/tolerance.band.md)
    — functional tolerance/prediction bands with methods: fpca,
    conformal, scb (Degras), exponential, elastic
- **Streaming depth**:
  - [`streaming.depth()`](https://sipemu.github.io/fdars-r/reference/streaming.depth.md)
    — fast depth computation against pre-sorted reference data
- **Soft-DTW distance**:
  - [`metric.softDTW()`](https://sipemu.github.io/fdars-r/reference/metric.softDTW.md)
    — differentiable Dynamic Time Warping
  - [`softdtw.barycenter()`](https://sipemu.github.io/fdars-r/reference/softdtw.barycenter.md)
    — Soft-DTW barycenter computation
- **Function-on-scalar regression**:
  - [`fosr()`](https://sipemu.github.io/fdars-r/reference/fosr.md) —
    penalized function-on-scalar regression
  - [`fosr.fpc()`](https://sipemu.github.io/fdars-r/reference/fosr.fpc.md)
    — FPC-based function-on-scalar regression
- **Functional ANOVA**:
  - [`fanova()`](https://sipemu.github.io/fdars-r/reference/fanova.md) —
    one-way functional ANOVA with permutation testing
- **Functional classification**:
  - [`fclassif()`](https://sipemu.github.io/fdars-r/reference/fclassif.md)
    — supervised classification (LDA, QDA, kNN, kernel, DD-plot)
  - [`fclassif.cv()`](https://sipemu.github.io/fdars-r/reference/fclassif.cv.md)
    — cross-validated classification error
- **Functional mixed models**:
  - [`fmm()`](https://sipemu.github.io/fdars-r/reference/fmm.md) —
    functional mixed models with subject-level random effects
  - [`fmm.predict()`](https://sipemu.github.io/fdars-r/reference/fmm.predict.md)
    — prediction from fitted FMM
  - [`fmm.test.fixed()`](https://sipemu.github.io/fdars-r/reference/fmm.test.fixed.md)
    — permutation test for fixed effects
- **GMM clustering**:
  - [`cluster.gmm()`](https://sipemu.github.io/fdars-r/reference/cluster.gmm.md)
    — Gaussian Mixture Model clustering with automatic K selection
  - [`cluster.init()`](https://sipemu.github.io/fdars-r/reference/cluster.init.md)
    — k-means++ initialization
- **Andrews transformation**:
  - [`andrews_transform()`](https://sipemu.github.io/fdars-r/reference/andrews_transform.md)
    — Fourier expansion of multivariate data to functional curves
  - [`andrews_loadings()`](https://sipemu.github.io/fdars-r/reference/andrews_loadings.md)
    — project FPCA eigenfunctions back to original variable space
- **Seasonal analysis** (comprehensive module):
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
- **Equivalence testing**:
  - [`fequiv.test()`](https://sipemu.github.io/fdars-r/reference/fequiv.test.md)
    — functional equivalence test
- **Functional linear model**:
  - [`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md),
    [`fregre.lm.cv()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.cv.md)
    — FPC-based functional linear regression
  - [`functional.logistic()`](https://sipemu.github.io/fdars-r/reference/functional.logistic.md)
    — functional logistic regression
- **Simulation toolbox**:
  - [`simFunData()`](https://sipemu.github.io/fdars-r/reference/simFunData.md),
    [`simMultiFunData()`](https://sipemu.github.io/fdars-r/reference/simMultiFunData.md)
    — simulate functional data via Karhunen-Loève
  - [`eFun()`](https://sipemu.github.io/fdars-r/reference/eFun.md),
    [`eVal()`](https://sipemu.github.io/fdars-r/reference/eVal.md) —
    eigenfunction/eigenvalue generators
  - [`addError()`](https://sipemu.github.io/fdars-r/reference/addError.md)
    — add measurement noise
  - [`irregFdata()`](https://sipemu.github.io/fdars-r/reference/irregFdata.md),
    [`sparsify()`](https://sipemu.github.io/fdars-r/reference/sparsify.md)
    — irregular/sparse functional data
- **Irregular functional data**:
  - [`irregFdata()`](https://sipemu.github.io/fdars-r/reference/irregFdata.md)
    — irregular functional data objects
  - [`as.fdata()`](https://sipemu.github.io/fdars-r/reference/as.fdata.irregFdata.md)
    — convert irregular to regular grid
  - [`mean.irregFdata()`](https://sipemu.github.io/fdars-r/reference/mean.irregFdata.md),
    [`summary.irregFdata()`](https://sipemu.github.io/fdars-r/reference/summary.irregFdata.md),
    [`plot.irregFdata()`](https://sipemu.github.io/fdars-r/reference/plot.irregFdata.md)
- **Unified dispatch interfaces**:
  - [`depth()`](https://sipemu.github.io/fdars-r/reference/depth.md) —
    dispatches to FM, mode, RP, RT, BD, MBD, MEI, FSD, KFSD, RPD
  - [`metric()`](https://sipemu.github.io/fdars-r/reference/metric.md) —
    dispatches to lp, hausdorff, dtw, softdtw, kl, elastic, amplitude,
    phase, pca, deriv, basis, fourier, hshift
- **Plot methods** for all new object types:
  - [`plot.elastic.align()`](https://sipemu.github.io/fdars-r/reference/plot.elastic.align.md),
    [`plot.karcher.mean()`](https://sipemu.github.io/fdars-r/reference/plot.karcher.mean.md),
    [`plot.alignment.quality()`](https://sipemu.github.io/fdars-r/reference/plot.alignment.quality.md)
  - [`plot.landmark.register()`](https://sipemu.github.io/fdars-r/reference/plot.landmark.register.md),
    [`plot.tsrvf()`](https://sipemu.github.io/fdars-r/reference/plot.tsrvf.md)
  - [`plot.tolerance.band()`](https://sipemu.github.io/fdars-r/reference/plot.tolerance.band.md)
  - [`plot.fclassif()`](https://sipemu.github.io/fdars-r/reference/plot.fclassif.md),
    [`plot.fmm()`](https://sipemu.github.io/fdars-r/reference/plot.fmm.md),
    [`plot.fosr()`](https://sipemu.github.io/fdars-r/reference/plot.fosr.md),
    [`plot.fanova()`](https://sipemu.github.io/fdars-r/reference/plot.fanova.md)
  - [`plot.cluster.gmm()`](https://sipemu.github.io/fdars-r/reference/plot.cluster.gmm.md)
- **Vignettes/articles** (17 articles covering all major features):
  - Elastic alignment, landmark registration, TSRVF, alignment
    comparison
  - Tolerance bands, streaming depth, distance metrics
  - Functional regression, classification, mixed models, GMM clustering
  - Andrews transformation, seasonal analysis
  - Example analyses: Berkeley growth study, wine dataset (4 articles)
- **Pre-built binary packages** for macOS (.tgz) and Windows (.zip)
  attached to GitHub Releases — install without Rust toolchain

### Changed

- Upgraded Rust backend (fdars-core) to v0.7.2
- All plots use ggplot2; respects `theme_set()`
- Renamed functions to avoid S3 dispatch conflicts:
  - `fdata2basis.cv()` →
    [`fdata2basis_cv()`](https://sipemu.github.io/fdars-r/reference/fdata2basis_cv.md)
  - `fdata2basis.2d()` →
    [`fdata2basis_2d()`](https://sipemu.github.io/fdars-r/reference/fdata2basis_2d.md)
  - `basis2fdata.2d()` →
    [`basis2fdata_2d()`](https://sipemu.github.io/fdars-r/reference/basis2fdata_2d.md)
  - [`norm.fdata()`](https://sipemu.github.io/fdars-r/reference/norm.md)
    → [`norm()`](https://sipemu.github.io/fdars-r/reference/norm.md)
  - `cov.*()` → `kernel.*()` (covariance kernel functions)
- Internal Rust wrapper functions no longer generate man pages (CRAN
  compliance)
- Portable Makefiles: removed GNU extensions (`export`, `:=`)
- Vendored Rust crates cleaned of hidden files and long paths

### Fixed

- `fdata2basis_2d` example error caused by extendr wrapper name
  collision
- `mean(fd)` now returns fdata object (was returning matrix)
- TSRVF inverse initial values
- Regression article cross-references and See Also sections

## \[0.4.0\] - 2024-12-13

### Added

- `id` and `metadata` slots in fdata objects for storing curve
  identifiers and associated data
- Outlier plot labeling: `plot(outliergram, label = "id")` or
  `label = "column_name"`
- [`magnitudeshape()`](https://sipemu.github.io/fdars-r/reference/magnitudeshape.md)
  labeling support (renamed from MS.plot)
- [`fregre.np.multi()`](https://sipemu.github.io/fdars-r/reference/fregre.np.multi.md)
  for regression with multiple functional predictors
- Enhanced
  [`plot.fdata()`](https://sipemu.github.io/fdars-r/reference/plot.fdata.md)
  with group coloring, mean curves, and confidence intervals
- [`group.distance()`](https://sipemu.github.io/fdars-r/reference/group.distance.md)
  for measuring distances between groups of curves
- [`group.test()`](https://sipemu.github.io/fdars-r/reference/group.test.md)
  permutation test for group differences
- [`outliergram()`](https://sipemu.github.io/fdars-r/reference/outliergram.md)
  visualization (MEI vs MBD plot)
- [`plot.fdata2pc()`](https://sipemu.github.io/fdars-r/reference/plot.fdata2pc.md)
  for FPCA visualization (components, variance, scores)

### Changed

- Auto-reduce alpha when `show.mean = TRUE` in
  [`plot.fdata()`](https://sipemu.github.io/fdars-r/reference/plot.fdata.md)
- Renamed `fdata.deriv()` to
  [`deriv()`](https://sipemu.github.io/fdars-r/reference/deriv.md) for
  consistency

### Fixed

- [`plot.group.distance()`](https://sipemu.github.io/fdars-r/reference/plot.group.distance.md)
  error handling
- Release workflow now generates documentation before building
- Missing `%||%` operator definition

## \[0.3.0\] - 2024-12-11

### Added

- Covariance kernel functions:
  [`kernel.exponential()`](https://sipemu.github.io/fdars-r/reference/kernel.exponential.md),
  [`kernel.matern()`](https://sipemu.github.io/fdars-r/reference/kernel.matern.md),
  [`kernel.gaussian()`](https://sipemu.github.io/fdars-r/reference/kernel.gaussian.md),
  etc.
- [`make.gaussian.process()`](https://sipemu.github.io/fdars-r/reference/make.gaussian.process.md)
  for simulating Gaussian process realizations
- 2D functional data support for most functions
- Unified API:
  [`depth()`](https://sipemu.github.io/fdars-r/reference/depth.md),
  [`median()`](https://sipemu.github.io/fdars-r/reference/median.md),
  [`trimmed()`](https://sipemu.github.io/fdars-r/reference/trimmed.md),
  [`trimvar()`](https://sipemu.github.io/fdars-r/reference/trimvar.md)
  with method parameter
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
  - [`fdata2basis()`](https://sipemu.github.io/fdars-r/reference/fdata2basis.md),
    [`basis2fdata()`](https://sipemu.github.io/fdars-r/reference/basis2fdata.md)
    — convert between data and basis coefficients
  - [`basis.gcv()`](https://sipemu.github.io/fdars-r/reference/basis.gcv.md),
    [`basis.aic()`](https://sipemu.github.io/fdars-r/reference/basis.aic.md),
    [`basis.bic()`](https://sipemu.github.io/fdars-r/reference/basis.bic.md)
    — goodness-of-fit metrics
  - [`fdata2basis_cv()`](https://sipemu.github.io/fdars-r/reference/fdata2basis_cv.md)
    — cross-validation for optimal nbasis selection
  - [`pspline()`](https://sipemu.github.io/fdars-r/reference/pspline.md)
    — P-spline smoothing with automatic lambda selection
  - [`fdata2basis_2d()`](https://sipemu.github.io/fdars-r/reference/fdata2basis_2d.md),
    [`basis2fdata_2d()`](https://sipemu.github.io/fdars-r/reference/basis2fdata_2d.md)
    — 2D tensor product basis support
  - [`pspline.2d()`](https://sipemu.github.io/fdars-r/reference/pspline.2d.md)
    — 2D P-spline with anisotropic penalties
  - [`select.basis.auto()`](https://sipemu.github.io/fdars-r/reference/select.basis.auto.md)
    — automatic basis selection per curve
  - [`normalize()`](https://sipemu.github.io/fdars-r/reference/normalize.md)
    — scale curves to unit Lp norm

### Changed

- Cleaned up API: removed backward compatibility shims
- Renamed functions for consistency (e.g., `fdata.mean` → `mean.fdata`)
- All plots now use ggplot2 instead of base R graphics

## \[0.1.0\] - 2024-12-10

### Initial Release

- **Core functional data structure**:
  [`fdata()`](https://sipemu.github.io/fdars-r/reference/fdata.md) for
  1D and 2D functional data
- **Depth functions**: FM, mode, RP, RT, FSD, KFSD, RPD (all with Rust
  backend)
- **Statistics**: mean, variance, standard deviation, covariance
- **Distance metrics**: Lp, Hausdorff, DTW, KL divergence
- **Semimetrics**: basis projection, Fourier, horizontal shift, PCA,
  derivative
- **Regression**:
  [`fregre.pc()`](https://sipemu.github.io/fdars-r/reference/fregre.pc.md),
  [`fregre.basis()`](https://sipemu.github.io/fdars-r/reference/fregre.basis.md),
  [`fregre.np()`](https://sipemu.github.io/fdars-r/reference/fregre.np.md)
  with CV variants
- **Outlier detection**:
  [`outliers.depth.pond()`](https://sipemu.github.io/fdars-r/reference/outliers.depth.pond.md),
  [`outliers.depth.trim()`](https://sipemu.github.io/fdars-r/reference/outliers.depth.trim.md),
  [`outliers.lrt()`](https://sipemu.github.io/fdars-r/reference/outliers.lrt.md)
- **Smoothing**: Nadaraya-Watson, local linear, local polynomial, k-NN
- **Clustering**: functional k-means
- **Hypothesis testing**:
  [`flm.test()`](https://sipemu.github.io/fdars-r/reference/flm.test.md),
  [`fmean.test.fdata()`](https://sipemu.github.io/fdars-r/reference/fmean.test.fdata.md)
- **Utilities**: Simpson integration, inner product, derivatives
