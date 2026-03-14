# Package index

## Functional Data Objects

- [`fdata()`](https://sipemu.github.io/fdars-r/reference/fdata.md) :
  Create a functional data object
- [`fdata.cen()`](https://sipemu.github.io/fdars-r/reference/fdata.cen.md)
  : Center functional data
- [`deriv()`](https://sipemu.github.io/fdars-r/reference/deriv.md) :
  Compute functional derivative
- [`fdata.gradient()`](https://sipemu.github.io/fdars-r/reference/fdata.gradient.md)
  : High-Accuracy Gradient for Functional Data
- [`mean(`*`<fdata>`*`)`](https://sipemu.github.io/fdars-r/reference/mean.fdata.md)
  : Compute functional mean
- [`median()`](https://sipemu.github.io/fdars-r/reference/median.md) :
  Compute Functional Median
- [`trimmed()`](https://sipemu.github.io/fdars-r/reference/trimmed.md) :
  Compute Functional Trimmed Mean
- [`trimvar()`](https://sipemu.github.io/fdars-r/reference/trimvar.md) :
  Compute Functional Trimmed Variance
- [`var()`](https://sipemu.github.io/fdars-r/reference/var.md) :
  Functional Variance
- [`sd()`](https://sipemu.github.io/fdars-r/reference/sd.md) :
  Functional Standard Deviation
- [`normalize()`](https://sipemu.github.io/fdars-r/reference/normalize.md)
  : Normalize functional data
- [`standardize()`](https://sipemu.github.io/fdars-r/reference/standardize.md)
  : Standardize functional data (z-score normalization)
- [`scale_minmax()`](https://sipemu.github.io/fdars-r/reference/scale_minmax.md)
  : Min-Max scaling for functional data
- [`gmed()`](https://sipemu.github.io/fdars-r/reference/gmed.md) :
  Geometric Median of Functional Data
- [`inprod.fdata()`](https://sipemu.github.io/fdars-r/reference/inprod.fdata.md)
  : Inner Product of Functional Data
- [`int.simpson()`](https://sipemu.github.io/fdars-r/reference/int.simpson.md)
  : Utility Functions for Functional Data Analysis
- [`localavg.fdata()`](https://sipemu.github.io/fdars-r/reference/localavg.fdata.md)
  : Local Averages Feature Extraction
- [`fdata.bootstrap()`](https://sipemu.github.io/fdars-r/reference/fdata.bootstrap.md)
  : Bootstrap Functional Data
- [`fdata.bootstrap.ci()`](https://sipemu.github.io/fdars-r/reference/fdata.bootstrap.ci.md)
  : Bootstrap Confidence Intervals for Functional Statistics
- [`df_to_fdata2d()`](https://sipemu.github.io/fdars-r/reference/df_to_fdata2d.md)
  : Convert DataFrame to 2D functional data

## Basis Representation

- [`fdata2basis()`](https://sipemu.github.io/fdars-r/reference/fdata2basis.md)
  : Convert Functional Data to Basis Coefficients
- [`fdata2basis_2d()`](https://sipemu.github.io/fdars-r/reference/fdata2basis_2d.md)
  : Convert 2D Functional Data to Tensor Product Basis Coefficients
- [`fdata2basis_cv()`](https://sipemu.github.io/fdars-r/reference/fdata2basis_cv.md)
  : Cross-Validation for Basis Function Number Selection
- [`basis2fdata()`](https://sipemu.github.io/fdars-r/reference/basis2fdata.md)
  : Basis Representation Functions for Functional Data
- [`basis2fdata_2d()`](https://sipemu.github.io/fdars-r/reference/basis2fdata_2d.md)
  : Reconstruct 2D Functional Data from Tensor Product Basis
  Coefficients
- [`fdata2fd()`](https://sipemu.github.io/fdars-r/reference/fdata2fd.md)
  : Convert Functional Data to fd class
- [`fdata2pc()`](https://sipemu.github.io/fdars-r/reference/fdata2pc.md)
  : Convert Functional Data to Principal Component Scores
- [`fdata2pls()`](https://sipemu.github.io/fdars-r/reference/fdata2pls.md)
  : Convert Functional Data to PLS Scores
- [`basis.aic()`](https://sipemu.github.io/fdars-r/reference/basis.aic.md)
  : AIC for Basis Representation
- [`basis.bic()`](https://sipemu.github.io/fdars-r/reference/basis.bic.md)
  : BIC for Basis Representation
- [`basis.gcv()`](https://sipemu.github.io/fdars-r/reference/basis.gcv.md)
  : GCV Score for Basis Representation
- [`select.basis.auto()`](https://sipemu.github.io/fdars-r/reference/select.basis.auto.md)
  : Automatic Per-Curve Basis Type and Number Selection
- [`pspline()`](https://sipemu.github.io/fdars-r/reference/pspline.md) :
  P-spline Smoothing for Functional Data
- [`pspline.2d()`](https://sipemu.github.io/fdars-r/reference/pspline.2d.md)
  : P-spline Smoothing for 2D Functional Data

## Andrews Transformation

- [`andrews_transform()`](https://sipemu.github.io/fdars-r/reference/andrews_transform.md)
  : Andrews Transformation
- [`andrews_loadings()`](https://sipemu.github.io/fdars-r/reference/andrews_loadings.md)
  : Andrews Loadings: Project FPCA Eigenfunctions to Original Variables

## Elastic Alignment

- [`srsf.transform()`](https://sipemu.github.io/fdars-r/reference/srsf.transform.md)
  : Elastic Alignment for Functional Data
- [`srsf.inverse()`](https://sipemu.github.io/fdars-r/reference/srsf.inverse.md)
  : Inverse SRSF Transform
- [`elastic.align()`](https://sipemu.github.io/fdars-r/reference/elastic.align.md)
  : Elastic Curve Alignment
- [`elastic.distance()`](https://sipemu.github.io/fdars-r/reference/elastic.distance.md)
  : Elastic Distance Matrix
- [`metric.elastic()`](https://sipemu.github.io/fdars-r/reference/metric.elastic.md)
  : Elastic Distance (Metric Dispatcher Alias)
- [`karcher.mean()`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md)
  : Karcher Mean in Elastic Metric
- [`periodic.rotate()`](https://sipemu.github.io/fdars-r/reference/periodic.rotate.md)
  : Periodic Rotation for Functional Data
- [`alignment.quality()`](https://sipemu.github.io/fdars-r/reference/alignment.quality.md)
  : Alignment Quality Diagnostics
- [`elastic.decomposition()`](https://sipemu.github.io/fdars-r/reference/elastic.decomposition.md)
  : Elastic Phase-Amplitude Decomposition
- [`amplitude.distance()`](https://sipemu.github.io/fdars-r/reference/amplitude.distance.md)
  : Amplitude Distance Matrix
- [`phase.distance()`](https://sipemu.github.io/fdars-r/reference/phase.distance.md)
  : Phase Distance Matrix
- [`elastic.align.constrained()`](https://sipemu.github.io/fdars-r/reference/elastic.align.constrained.md)
  : Landmark-Constrained Elastic Alignment
- [`alignment.pairwise.consistency()`](https://sipemu.github.io/fdars-r/reference/alignment.pairwise.consistency.md)
  : Pairwise Alignment Consistency
- [`elastic.pair()`](https://sipemu.github.io/fdars-r/reference/elastic.pair.md)
  : Elastic Pairwise Alignment
- [`srsf.reparameterize()`](https://sipemu.github.io/fdars-r/reference/srsf.reparameterize.md)
  : Apply Warping Function to a Curve
- [`warp.complexity()`](https://sipemu.github.io/fdars-r/reference/warp.complexity.md)
  : Warping Complexity
- [`warp.compose()`](https://sipemu.github.io/fdars-r/reference/warp.compose.md)
  : Compose Two Warping Functions
- [`warp.smoothness()`](https://sipemu.github.io/fdars-r/reference/warp.smoothness.md)
  : Warping Smoothness
- [`plot(`*`<elastic.align>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.elastic.align.md)
  : Plot Elastic Alignment Results
- [`plot(`*`<karcher.mean>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.karcher.mean.md)
  : Plot Karcher Mean Results
- [`plot(`*`<alignment.quality>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.alignment.quality.md)
  : Plot Alignment Quality Diagnostics

## Landmark Registration

- [`detect.landmarks()`](https://sipemu.github.io/fdars-r/reference/detect.landmarks.md)
  : Landmark Registration for Functional Data
- [`landmark.register()`](https://sipemu.github.io/fdars-r/reference/landmark.register.md)
  : Landmark Registration
- [`plot(`*`<landmark.register>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.landmark.register.md)
  : Plot Landmark Registration Results

## TSRVF

- [`tsrvf.transform()`](https://sipemu.github.io/fdars-r/reference/tsrvf.transform.md)
  : TSRVF: Transported Square-Root Velocity Function
- [`tsrvf.from.alignment()`](https://sipemu.github.io/fdars-r/reference/tsrvf.from.alignment.md)
  : TSRVF from Pre-computed Alignment
- [`tsrvf.inverse()`](https://sipemu.github.io/fdars-r/reference/tsrvf.inverse.md)
  : Inverse TSRVF Transform
- [`plot(`*`<tsrvf>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.tsrvf.md)
  : Plot TSRVF Results

## Tolerance Bands

- [`tolerance.band()`](https://sipemu.github.io/fdars-r/reference/tolerance.band.md)
  : Tolerance Bands for Functional Data
- [`plot(`*`<tolerance.band>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.tolerance.band.md)
  : Plot Tolerance Band

## Depth Functions

- [`depth()`](https://sipemu.github.io/fdars-r/reference/depth.md) :
  Depth Functions for Functional Data
- [`depth.BD()`](https://sipemu.github.io/fdars-r/reference/depth.BD.md)
  : Band Depth
- [`depth.FM()`](https://sipemu.github.io/fdars-r/reference/depth.FM.md)
  : Fraiman-Muniz Depth
- [`depth.FSD()`](https://sipemu.github.io/fdars-r/reference/depth.FSD.md)
  : Functional Spatial Depth
- [`depth.KFSD()`](https://sipemu.github.io/fdars-r/reference/depth.KFSD.md)
  : Kernel Functional Spatial Depth
- [`depth.MBD()`](https://sipemu.github.io/fdars-r/reference/depth.MBD.md)
  : Modified Band Depth
- [`depth.MEI()`](https://sipemu.github.io/fdars-r/reference/depth.MEI.md)
  : Modified Epigraph Index
- [`depth.mode()`](https://sipemu.github.io/fdars-r/reference/depth.mode.md)
  : Modal Depth
- [`depth.RP()`](https://sipemu.github.io/fdars-r/reference/depth.RP.md)
  : Random Projection Depth
- [`rp.stat()`](https://sipemu.github.io/fdars-r/reference/rp.stat.md) :
  Random Projection Statistic
- [`depth.RPD()`](https://sipemu.github.io/fdars-r/reference/depth.RPD.md)
  : Random Projection Depth with Derivatives
- [`depth.RT()`](https://sipemu.github.io/fdars-r/reference/depth.RT.md)
  : Random Tukey Depth
- [`streaming.depth()`](https://sipemu.github.io/fdars-r/reference/streaming.depth.md)
  : Streaming Depth for Functional Data
- [`depth.streaming()`](https://sipemu.github.io/fdars-r/reference/depth.streaming.md)
  : Streaming Depth (Alias)

## Distance & Metrics

- [`metric()`](https://sipemu.github.io/fdars-r/reference/metric.md) :
  Distance Metrics for Functional Data
- [`metric.DTW()`](https://sipemu.github.io/fdars-r/reference/metric.DTW.md)
  : Dynamic Time Warping for Functional Data
- [`metric.hausdorff()`](https://sipemu.github.io/fdars-r/reference/metric.hausdorff.md)
  : Hausdorff Metric for Functional Data
- [`metric.kl()`](https://sipemu.github.io/fdars-r/reference/metric.kl.md)
  : Kullback-Leibler Divergence Metric for Functional Data
- [`metric.lp()`](https://sipemu.github.io/fdars-r/reference/metric.lp.md)
  : Lp Metric for Functional Data
- [`metric.softDTW()`](https://sipemu.github.io/fdars-r/reference/metric.softDTW.md)
  : Soft Dynamic Time Warping Distance
- [`softdtw.barycenter()`](https://sipemu.github.io/fdars-r/reference/softdtw.barycenter.md)
  : Soft-DTW Barycenter
- [`norm()`](https://sipemu.github.io/fdars-r/reference/norm.md) :
  Compute Lp Norm of Functional Data
- [`semimetric.basis()`](https://sipemu.github.io/fdars-r/reference/semimetric.basis.md)
  : Semi-metric based on Basis Expansion
- [`semimetric.deriv()`](https://sipemu.github.io/fdars-r/reference/semimetric.deriv.md)
  : Semi-metric based on Derivatives
- [`semimetric.fourier()`](https://sipemu.github.io/fdars-r/reference/semimetric.fourier.md)
  : Semi-metric based on Fourier Coefficients (FFT)
- [`semimetric.hshift()`](https://sipemu.github.io/fdars-r/reference/semimetric.hshift.md)
  : Semi-metric based on Horizontal Shift (Time Warping)
- [`semimetric.pca()`](https://sipemu.github.io/fdars-r/reference/semimetric.pca.md)
  : Semi-metric based on Principal Components
- [`group.distance()`](https://sipemu.github.io/fdars-r/reference/group.distance.md)
  : Compute Distance/Similarity Between Groups of Functional Data

## Clustering

- [`cluster.fcm()`](https://sipemu.github.io/fdars-r/reference/cluster.fcm.md)
  : Fuzzy C-Means Clustering for Functional Data
- [`cluster.gmm()`](https://sipemu.github.io/fdars-r/reference/cluster.gmm.md)
  : Gaussian Mixture Model Clustering for Functional Data
- [`gmm.em()`](https://sipemu.github.io/fdars-r/reference/gmm.em.md) :
  Gaussian Mixture Model EM on Feature Matrix
- [`cluster.init()`](https://sipemu.github.io/fdars-r/reference/cluster.init.md)
  : K-Means++ Center Initialization
- [`cluster.kmeans()`](https://sipemu.github.io/fdars-r/reference/cluster.kmeans.md)
  : Clustering Functions for Functional Data
- [`cluster.optim()`](https://sipemu.github.io/fdars-r/reference/cluster.optim.md)
  : Optimal Number of Clusters for Functional K-Means

## Outlier Detection

- [`outliergram()`](https://sipemu.github.io/fdars-r/reference/outliergram.md)
  : Outliergram for Functional Data
- [`outliers.boxplot()`](https://sipemu.github.io/fdars-r/reference/outliers.boxplot.md)
  : Outlier Detection using Functional Boxplot
- [`outliers.depth.pond()`](https://sipemu.github.io/fdars-r/reference/outliers.depth.pond.md)
  : Outlier Detection for Functional Data
- [`outliers.depth.trim()`](https://sipemu.github.io/fdars-r/reference/outliers.depth.trim.md)
  : Outlier Detection using Trimmed Depth
- [`outliers.lrt()`](https://sipemu.github.io/fdars-r/reference/outliers.lrt.md)
  : LRT-based Outlier Detection for Functional Data
- [`outliers.thres.lrt()`](https://sipemu.github.io/fdars-r/reference/outliers.thres.lrt.md)
  : LRT Outlier Detection Threshold
- [`outliers.lrt.dist()`](https://sipemu.github.io/fdars-r/reference/outliers.lrt.dist.md)
  : LRT Outlier Threshold with Bootstrap Distribution
- [`magnitudeshape()`](https://sipemu.github.io/fdars-r/reference/magnitudeshape.md)
  : Magnitude-Shape Outlier Detection for Functional Data

## Regression

- [`fregre.basis()`](https://sipemu.github.io/fdars-r/reference/fregre.basis.md)
  : Functional Basis Regression
- [`fregre.basis.cv()`](https://sipemu.github.io/fdars-r/reference/fregre.basis.cv.md)
  : Cross-Validation for Functional Basis Regression
- [`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md)
  : Functional Linear Model (FPC-based)
- [`fregre.lm.cv()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.cv.md)
  : Cross-Validation for FPC Component Selection (fregre.lm)
- [`fregre.bootstrap.ci()`](https://sipemu.github.io/fdars-r/reference/fregre.bootstrap.ci.md)
  : Bootstrap Confidence Intervals for Functional Coefficient
- [`fregre.np()`](https://sipemu.github.io/fdars-r/reference/fregre.np.md)
  : Nonparametric Functional Regression
- [`fregre.np.cv()`](https://sipemu.github.io/fdars-r/reference/fregre.np.cv.md)
  : Cross-Validation for Nonparametric Functional Regression
- [`fregre.np.mixed()`](https://sipemu.github.io/fdars-r/reference/fregre.np.mixed.md)
  : Nonparametric Functional Regression with Mixed Predictors
- [`fregre.np.multi()`](https://sipemu.github.io/fdars-r/reference/fregre.np.multi.md)
  : Nonparametric Regression with Multiple Functional Predictors
- [`fregre.pc()`](https://sipemu.github.io/fdars-r/reference/fregre.pc.md)
  : Functional Regression
- [`fregre.pc.cv()`](https://sipemu.github.io/fdars-r/reference/fregre.pc.cv.md)
  : Cross-Validation for Functional PC Regression
- [`functional.logistic()`](https://sipemu.github.io/fdars-r/reference/functional.logistic.md)
  : Functional Logistic Regression
- [`predict(`*`<fregre.logistic>`*`)`](https://sipemu.github.io/fdars-r/reference/predict.fregre.logistic.md)
  : Predict from Functional Logistic Model
- [`model.selection.ncomp()`](https://sipemu.github.io/fdars-r/reference/model.selection.ncomp.md)
  : Model Selection for Number of FPC Components
- [`optim.np()`](https://sipemu.github.io/fdars-r/reference/optim.np.md)
  : Optimize Bandwidth Using Cross-Validation
- [`flm.test()`](https://sipemu.github.io/fdars-r/reference/flm.test.md)
  : Statistical Tests for Functional Data
- [`pred.MAE()`](https://sipemu.github.io/fdars-r/reference/pred.MAE.md)
  : Mean Absolute Error
- [`pred.MSE()`](https://sipemu.github.io/fdars-r/reference/pred.MSE.md)
  : Mean Squared Error
- [`pred.R2()`](https://sipemu.github.io/fdars-r/reference/pred.R2.md) :
  R-Squared (Coefficient of Determination)
- [`pred.RMSE()`](https://sipemu.github.io/fdars-r/reference/pred.RMSE.md)
  : Root Mean Squared Error

## Function-on-Scalar Regression

- [`fosr()`](https://sipemu.github.io/fdars-r/reference/fosr.md) :
  Function-on-Scalar Regression
- [`fosr.fpc()`](https://sipemu.github.io/fdars-r/reference/fosr.fpc.md)
  : FPC-based Function-on-Scalar Regression
- [`fosr.2d()`](https://sipemu.github.io/fdars-r/reference/fosr.2d.md) :
  2D Function-on-Scalar Regression
- [`fanova()`](https://sipemu.github.io/fdars-r/reference/fanova.md) :
  Functional ANOVA

## Classification

- [`fclassif()`](https://sipemu.github.io/fdars-r/reference/fclassif.md)
  : Supervised Classification of Functional Data
- [`fclassif.cv()`](https://sipemu.github.io/fdars-r/reference/fclassif.cv.md)
  : Cross-Validated Functional Classification

## Cross-Validation

- [`cv.fdata()`](https://sipemu.github.io/fdars-r/reference/cv.fdata.md)
  : Unified K-Fold Cross-Validation for Functional Data
- [`plot(`*`<cv.fdata>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.cv.fdata.md)
  : Plot Method for cv.fdata
- [`print(`*`<cv.fdata>`*`)`](https://sipemu.github.io/fdars-r/reference/print.cv.fdata.md)
  : Print Method for cv.fdata

## Functional Mixed Models

- [`fmm()`](https://sipemu.github.io/fdars-r/reference/fmm.md) :
  Functional Mixed Models
- [`fmm.predict()`](https://sipemu.github.io/fdars-r/reference/fmm.predict.md)
  : Predict from Functional Mixed Model
- [`fmm.test.fixed()`](https://sipemu.github.io/fdars-r/reference/fmm.test.fixed.md)
  : Permutation Test for Fixed Effects in FMM

## Seasonal Analysis

- [`estimate.period()`](https://sipemu.github.io/fdars-r/reference/estimate.period.md)
  : Estimate Seasonal Period using FFT
- [`detect.period()`](https://sipemu.github.io/fdars-r/reference/detect.period.md)
  : Seasonal Analysis Functions for Functional Data
- [`detect.periods()`](https://sipemu.github.io/fdars-r/reference/detect.periods.md)
  : Detect Multiple Concurrent Periods
- [`detect.peaks()`](https://sipemu.github.io/fdars-r/reference/detect.peaks.md)
  : Detect Peaks in Functional Data
- [`autoperiod()`](https://sipemu.github.io/fdars-r/reference/autoperiod.md)
  : Autoperiod: Hybrid FFT + ACF Period Detection
- [`cfd.autoperiod()`](https://sipemu.github.io/fdars-r/reference/cfd.autoperiod.md)
  : CFDAutoperiod: Clustered Filtered Detrended Autoperiod
- [`sazed()`](https://sipemu.github.io/fdars-r/reference/sazed.md) :
  SAZED: Spectral-ACF Zero-crossing Ensemble Detection
- [`lomb.scargle()`](https://sipemu.github.io/fdars-r/reference/lomb.scargle.md)
  : Lomb-Scargle Periodogram
- [`matrix.profile()`](https://sipemu.github.io/fdars-r/reference/matrix.profile.md)
  : Matrix Profile for Motif Discovery and Period Detection
- [`stl.fd()`](https://sipemu.github.io/fdars-r/reference/stl.fd.md) :
  STL Decomposition: Seasonal and Trend decomposition using LOESS
- [`ssa.fd()`](https://sipemu.github.io/fdars-r/reference/ssa.fd.md) :
  Singular Spectrum Analysis (SSA) for Time Series Decomposition
- [`seasonal.strength()`](https://sipemu.github.io/fdars-r/reference/seasonal.strength.md)
  : Measure Seasonal Strength
- [`seasonal.strength.curve()`](https://sipemu.github.io/fdars-r/reference/seasonal.strength.curve.md)
  : Time-Varying Seasonal Strength
- [`detect.seasonality.changes()`](https://sipemu.github.io/fdars-r/reference/detect.seasonality.changes.md)
  : Detect Changes in Seasonality
- [`detect.seasonality.changes.auto()`](https://sipemu.github.io/fdars-r/reference/detect.seasonality.changes.auto.md)
  : Detect Seasonality Changes with Automatic Threshold
- [`detect_amplitude_modulation()`](https://sipemu.github.io/fdars-r/reference/detect_amplitude_modulation.md)
  : Detect Amplitude Modulation in Seasonal Time Series
- [`instantaneous.period()`](https://sipemu.github.io/fdars-r/reference/instantaneous.period.md)
  : Estimate Instantaneous Period
- [`analyze.peak.timing()`](https://sipemu.github.io/fdars-r/reference/analyze.peak.timing.md)
  : Analyze Peak Timing Variability
- [`classify.seasonality()`](https://sipemu.github.io/fdars-r/reference/classify.seasonality.md)
  : Classify Seasonality Type
- [`detrend()`](https://sipemu.github.io/fdars-r/reference/detrend.md) :
  Remove Trend from Functional Data
- [`decompose()`](https://sipemu.github.io/fdars-r/reference/decompose.md)
  : Seasonal-Trend Decomposition

## Penalized Basis Smoothing

- [`smooth.basis.fd()`](https://sipemu.github.io/fdars-r/reference/smooth.basis.fd.md)
  : Penalized Basis Smoothing
- [`smooth.basis.gcv()`](https://sipemu.github.io/fdars-r/reference/smooth.basis.gcv.md)
  : Penalized Basis Smoothing with GCV-Optimal Lambda

## Elastic FPCA

- [`vert.fpca()`](https://sipemu.github.io/fdars-r/reference/vert.fpca.md)
  : Elastic FPCA
- [`horiz.fpca()`](https://sipemu.github.io/fdars-r/reference/horiz.fpca.md)
  : Horizontal (Phase) FPCA
- [`joint.fpca()`](https://sipemu.github.io/fdars-r/reference/joint.fpca.md)
  : Joint (Amplitude + Phase) FPCA

## Elastic Regression

- [`elastic.regression()`](https://sipemu.github.io/fdars-r/reference/elastic.regression.md)
  : Elastic Regression
- [`elastic.logistic()`](https://sipemu.github.io/fdars-r/reference/elastic.logistic.md)
  : Elastic Logistic Classification
- [`elastic.pcr()`](https://sipemu.github.io/fdars-r/reference/elastic.pcr.md)
  : Elastic Principal Component Regression
- [`predict(`*`<elastic.regression>`*`)`](https://sipemu.github.io/fdars-r/reference/predict.elastic.regression.md)
  : Predict from Elastic Regression
- [`predict(`*`<elastic.logistic>`*`)`](https://sipemu.github.io/fdars-r/reference/predict.elastic.logistic.md)
  : Predict from Elastic Logistic Classification

## Scalar-on-Shape Regression

- [`scalar.on.shape()`](https://sipemu.github.io/fdars-r/reference/scalar.on.shape.md)
  : Scalar-on-Shape Regression
- [`predict(`*`<scalar.on.shape>`*`)`](https://sipemu.github.io/fdars-r/reference/predict.scalar.on.shape.md)
  : Predict from a Scalar-on-Shape Model
- [`print(`*`<scalar.on.shape>`*`)`](https://sipemu.github.io/fdars-r/reference/print.scalar.on.shape.md)
  : Print Scalar-on-Shape Model

## Robust Regression

- [`fregre.l1()`](https://sipemu.github.io/fdars-r/reference/fregre.l1.md)
  : L1 (Least Absolute Deviation) Functional Regression
- [`fregre.huber()`](https://sipemu.github.io/fdars-r/reference/fregre.huber.md)
  : Huber M-Estimation Functional Regression
- [`predict(`*`<fregre.robust>`*`)`](https://sipemu.github.io/fdars-r/reference/predict.fregre.robust.md)
  : Predict from Robust Functional Regression

## Statistical Process Monitoring

- [`spm`](https://sipemu.github.io/fdars-r/reference/spm.md) :
  Statistical Process Monitoring for Functional Data
- [`spm.phase1()`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md)
  : Build Univariate SPM Control Chart (Phase I)
- [`spm.monitor()`](https://sipemu.github.io/fdars-r/reference/spm.monitor.md)
  : Monitor New Functional Data (Phase II)
- [`spm.ewma()`](https://sipemu.github.io/fdars-r/reference/spm.ewma.md)
  : EWMA-Based SPM Monitoring
- [`spm.contributions()`](https://sipemu.github.io/fdars-r/reference/spm.contributions.md)
  : SPM Contribution Diagnostics
- [`mfpca()`](https://sipemu.github.io/fdars-r/reference/mfpca.md) :
  Multivariate Functional Principal Component Analysis
- [`frcc.phase1()`](https://sipemu.github.io/fdars-r/reference/frcc.phase1.md)
  : Build Functional Regression Control Chart (Phase I)
- [`frcc.monitor()`](https://sipemu.github.io/fdars-r/reference/frcc.monitor.md)
  : Monitor New Data Against FRCC Chart (Phase II)

## Elastic Changepoint

- [`elastic.changepoint()`](https://sipemu.github.io/fdars-r/reference/elastic.changepoint.md)
  : Elastic Changepoint Detection

## Conformal Prediction

- [`conformal.fregre.lm()`](https://sipemu.github.io/fdars-r/reference/conformal.fregre.lm.md)
  : Conformal Prediction for Functional Linear Model
- [`conformal.fregre.np()`](https://sipemu.github.io/fdars-r/reference/conformal.fregre.np.md)
  : Conformal Prediction for Nonparametric Functional Regression
- [`conformal.elastic.regression()`](https://sipemu.github.io/fdars-r/reference/conformal.elastic.regression.md)
  : Conformal Prediction for Elastic Regression
- [`conformal.elastic.pcr()`](https://sipemu.github.io/fdars-r/reference/conformal.elastic.pcr.md)
  : Conformal Prediction for Elastic PCR
- [`conformal.elastic.logistic()`](https://sipemu.github.io/fdars-r/reference/conformal.elastic.logistic.md)
  : Conformal Prediction for Elastic Logistic Regression
- [`conformal.logistic()`](https://sipemu.github.io/fdars-r/reference/conformal.logistic.md)
  : Conformal Prediction for Logistic Regression
- [`conformal.classif()`](https://sipemu.github.io/fdars-r/reference/conformal.classif.md)
  : Conformal Classification
- [`cv.conformal.regression()`](https://sipemu.github.io/fdars-r/reference/cv.conformal.regression.md)
  : Cross-Conformal (CV+) Regression
- [`cv.conformal.classification()`](https://sipemu.github.io/fdars-r/reference/cv.conformal.classification.md)
  : Cross-Conformal (CV+) Classification
- [`jackknife.plus()`](https://sipemu.github.io/fdars-r/reference/jackknife.plus.md)
  : Jackknife+ Regression
- [`conformal.generic.regression()`](https://sipemu.github.io/fdars-r/reference/conformal.generic.regression.md)
  : Generic Conformal Regression
- [`conformal.generic.classification()`](https://sipemu.github.io/fdars-r/reference/conformal.generic.classification.md)
  : Generic Conformal Classification

## Explainability

- [`fregre.pdp()`](https://sipemu.github.io/fdars-r/reference/fregre.pdp.md)
  : Functional Partial Dependence Plot
- [`fregre.shap()`](https://sipemu.github.io/fdars-r/reference/fregre.shap.md)
  : FPC SHAP Values
- [`fregre.ale()`](https://sipemu.github.io/fdars-r/reference/fregre.ale.md)
  : Accumulated Local Effects
- [`fregre.lime()`](https://sipemu.github.io/fdars-r/reference/fregre.lime.md)
  : LIME Explanation
- [`fregre.anchor()`](https://sipemu.github.io/fdars-r/reference/fregre.anchor.md)
  : Anchor Explanation
- [`fregre.counterfactual()`](https://sipemu.github.io/fdars-r/reference/fregre.counterfactual.md)
  : Counterfactual Explanation
- [`fregre.saliency()`](https://sipemu.github.io/fdars-r/reference/fregre.saliency.md)
  : Functional Saliency Map
- [`fregre.importance()`](https://sipemu.github.io/fdars-r/reference/fregre.importance.md)
  : FPC Permutation Importance
- [`fregre.conditional.importance()`](https://sipemu.github.io/fdars-r/reference/fregre.conditional.importance.md)
  : Conditional Permutation Importance
- [`fregre.influence()`](https://sipemu.github.io/fdars-r/reference/fregre.influence.md)
  : Influence Diagnostics
- [`fregre.vif()`](https://sipemu.github.io/fdars-r/reference/fregre.vif.md)
  : Variance Inflation Factors
- [`fregre.dfbetas()`](https://sipemu.github.io/fdars-r/reference/fregre.dfbetas.md)
  : DFBETAS and DFFITS
- [`fregre.loo()`](https://sipemu.github.io/fdars-r/reference/fregre.loo.md)
  : LOO-CV and PRESS
- [`fregre.sobol()`](https://sipemu.github.io/fdars-r/reference/fregre.sobol.md)
  : Sobol Indices
- [`fregre.friedman()`](https://sipemu.github.io/fdars-r/reference/fregre.friedman.md)
  : Friedman H-Statistic
- [`fregre.conformal()`](https://sipemu.github.io/fdars-r/reference/fregre.conformal.md)
  : Conformal Prediction
- [`fregre.stability()`](https://sipemu.github.io/fdars-r/reference/fregre.stability.md)
  : Explanation Stability
- [`fregre.depth()`](https://sipemu.github.io/fdars-r/reference/fregre.depth.md)
  : Regression Depth
- [`fregre.domain()`](https://sipemu.github.io/fdars-r/reference/fregre.domain.md)
  : Domain Selection
- [`fregre.prediction.interval()`](https://sipemu.github.io/fdars-r/reference/fregre.prediction.interval.md)
  : Prediction Intervals
- [`fregre.prototype()`](https://sipemu.github.io/fdars-r/reference/fregre.prototype.md)
  : Prototype and Criticism
- [`fregre.beta.decomp()`](https://sipemu.github.io/fdars-r/reference/fregre.beta.decomp.md)
  : Beta Decomposition
- [`fregre.pointwise()`](https://sipemu.github.io/fdars-r/reference/fregre.pointwise.md)
  : Pointwise Importance
- [`fregre.significant.regions()`](https://sipemu.github.io/fdars-r/reference/fregre.significant.regions.md)
  : Significant Regions
- [`fregre.calibration()`](https://sipemu.github.io/fdars-r/reference/fregre.calibration.md)
  : Calibration Diagnostics (Logistic)
- [`fregre.ece()`](https://sipemu.github.io/fdars-r/reference/fregre.ece.md)
  : Expected Calibration Error (Logistic)
- [`elastic.attribution()`](https://sipemu.github.io/fdars-r/reference/elastic.attribution.md)
  : Elastic PCR Attribution

## Smoothing

- [`S.KNN()`](https://sipemu.github.io/fdars-r/reference/S.KNN.md) :
  K-Nearest Neighbors Smoother Matrix
- [`S.LCR()`](https://sipemu.github.io/fdars-r/reference/S.LCR.md) :
  Local Cubic Regression Smoother Matrix
- [`S.LLR()`](https://sipemu.github.io/fdars-r/reference/S.LLR.md) :
  Local Linear Regression Smoother Matrix
- [`S.LPR()`](https://sipemu.github.io/fdars-r/reference/S.LPR.md) :
  Local Polynomial Regression Smoother Matrix
- [`S.NW()`](https://sipemu.github.io/fdars-r/reference/S.NW.md) :
  Smoothing Functions for Functional Data
- [`CV.S()`](https://sipemu.github.io/fdars-r/reference/CV.S.md) :
  Cross-Validation for Smoother Selection
- [`GCV.S()`](https://sipemu.github.io/fdars-r/reference/GCV.S.md) :
  Generalized Cross-Validation for Smoother Selection
- [`h.default()`](https://sipemu.github.io/fdars-r/reference/h.default.md)
  : Default Bandwidth
- [`register.fd()`](https://sipemu.github.io/fdars-r/reference/register.fd.md)
  : Curve Registration (Alignment)

## Kernels (Smoothing)

- [`Kernel()`](https://sipemu.github.io/fdars-r/reference/Kernel.md) :
  Unified Symmetric Kernel Interface
- [`Kernel.asymmetric()`](https://sipemu.github.io/fdars-r/reference/Kernel.asymmetric.md)
  : Unified Asymmetric Kernel Interface
- [`Kernel.integrate()`](https://sipemu.github.io/fdars-r/reference/Kernel.integrate.md)
  : Unified Integrated Kernel Interface
- [`Ker.cos()`](https://sipemu.github.io/fdars-r/reference/Ker.cos.md) :
  Cosine Kernel
- [`Ker.epa()`](https://sipemu.github.io/fdars-r/reference/Ker.epa.md) :
  Epanechnikov Kernel
- [`Ker.norm()`](https://sipemu.github.io/fdars-r/reference/Ker.norm.md)
  : Kernel Functions
- [`Ker.quar()`](https://sipemu.github.io/fdars-r/reference/Ker.quar.md)
  : Quartic (Biweight) Kernel
- [`Ker.tri()`](https://sipemu.github.io/fdars-r/reference/Ker.tri.md) :
  Triweight Kernel
- [`Ker.unif()`](https://sipemu.github.io/fdars-r/reference/Ker.unif.md)
  : Uniform (Rectangular) Kernel
- [`AKer.cos()`](https://sipemu.github.io/fdars-r/reference/AKer.cos.md)
  : Asymmetric Cosine Kernel
- [`AKer.epa()`](https://sipemu.github.io/fdars-r/reference/AKer.epa.md)
  : Asymmetric Epanechnikov Kernel
- [`AKer.norm()`](https://sipemu.github.io/fdars-r/reference/AKer.norm.md)
  : Asymmetric Normal Kernel
- [`AKer.quar()`](https://sipemu.github.io/fdars-r/reference/AKer.quar.md)
  : Asymmetric Quartic Kernel
- [`AKer.tri()`](https://sipemu.github.io/fdars-r/reference/AKer.tri.md)
  : Asymmetric Triweight Kernel
- [`AKer.unif()`](https://sipemu.github.io/fdars-r/reference/AKer.unif.md)
  : Asymmetric Uniform Kernel
- [`IKer.cos()`](https://sipemu.github.io/fdars-r/reference/IKer.cos.md)
  : Integrated Cosine Kernel
- [`IKer.epa()`](https://sipemu.github.io/fdars-r/reference/IKer.epa.md)
  : Integrated Epanechnikov Kernel
- [`IKer.norm()`](https://sipemu.github.io/fdars-r/reference/IKer.norm.md)
  : Integrated Normal Kernel
- [`IKer.quar()`](https://sipemu.github.io/fdars-r/reference/IKer.quar.md)
  : Integrated Quartic Kernel
- [`IKer.tri()`](https://sipemu.github.io/fdars-r/reference/IKer.tri.md)
  : Integrated Triweight Kernel
- [`IKer.unif()`](https://sipemu.github.io/fdars-r/reference/IKer.unif.md)
  : Integrated Uniform Kernel

## Covariance Functions (GP)

- [`kernel.add()`](https://sipemu.github.io/fdars-r/reference/kernel.add.md)
  : Add Covariance Functions
- [`kernel.brownian()`](https://sipemu.github.io/fdars-r/reference/kernel.brownian.md)
  : Brownian Motion Covariance Function
- [`kernel.exponential()`](https://sipemu.github.io/fdars-r/reference/kernel.exponential.md)
  : Exponential Covariance Function
- [`kernel.gaussian()`](https://sipemu.github.io/fdars-r/reference/kernel.gaussian.md)
  : Gaussian (Squared Exponential) Covariance Function
- [`kernel.linear()`](https://sipemu.github.io/fdars-r/reference/kernel.linear.md)
  : Linear Covariance Function
- [`kernel.matern()`](https://sipemu.github.io/fdars-r/reference/kernel.matern.md)
  : Matern Covariance Function
- [`kernel.mult()`](https://sipemu.github.io/fdars-r/reference/kernel.mult.md)
  : Multiply Covariance Functions
- [`kernel.periodic()`](https://sipemu.github.io/fdars-r/reference/kernel.periodic.md)
  : Periodic Covariance Function
- [`kernel.polynomial()`](https://sipemu.github.io/fdars-r/reference/kernel.polynomial.md)
  : Polynomial Covariance Function
- [`kernel.whitenoise()`](https://sipemu.github.io/fdars-r/reference/kernel.whitenoise.md)
  : White Noise Covariance Function
- [`make.gaussian.process()`](https://sipemu.github.io/fdars-r/reference/make.gaussian.process.md)
  : Generate Gaussian Process Samples
- [`cov()`](https://sipemu.github.io/fdars-r/reference/cov.md) :
  Functional Covariance Function

## Simulation

- [`eFun()`](https://sipemu.github.io/fdars-r/reference/eFun.md) :
  Generate Eigenfunction Basis
- [`eVal()`](https://sipemu.github.io/fdars-r/reference/eVal.md) :
  Generate Eigenvalue Sequence
- [`simFunData()`](https://sipemu.github.io/fdars-r/reference/simFunData.md)
  : Simulate Functional Data via Karhunen-Loeve Expansion
- [`simMultiFunData()`](https://sipemu.github.io/fdars-r/reference/simMultiFunData.md)
  : Simulate Multivariate Functional Data
- [`addError()`](https://sipemu.github.io/fdars-r/reference/addError.md)
  : Add Measurement Error to Functional Data

## Irregular Functional Data

- [`irregFdata()`](https://sipemu.github.io/fdars-r/reference/irregFdata.md)
  : Create an Irregular Functional Data Object
- [`is.irregular()`](https://sipemu.github.io/fdars-r/reference/is.irregular.md)
  : Check if an Object is Irregular Functional Data
- [`sparsify()`](https://sipemu.github.io/fdars-r/reference/sparsify.md)
  : Convert Regular Functional Data to Irregular by Subsampling
- [`as.fdata()`](https://sipemu.github.io/fdars-r/reference/as.fdata.irregFdata.md)
  : Convert Irregular Functional Data to Regular Grid
- [`mean(`*`<irregFdata>`*`)`](https://sipemu.github.io/fdars-r/reference/mean.irregFdata.md)
  : Estimate Mean Function for Irregular Data
- [`summary(`*`<irregFdata>`*`)`](https://sipemu.github.io/fdars-r/reference/summary.irregFdata.md)
  : Summary method for irregFdata objects
- [`print(`*`<irregFdata>`*`)`](https://sipemu.github.io/fdars-r/reference/print.irregFdata.md)
  : Print method for irregFdata objects
- [`autoplot(`*`<irregFdata>`*`)`](https://sipemu.github.io/fdars-r/reference/autoplot.irregFdata.md)
  : Autoplot method for irregFdata objects
- [`plot(`*`<irregFdata>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.irregFdata.md)
  : Plot method for irregFdata objects
- [`` `[`( ``*`<irregFdata>`*`)`](https://sipemu.github.io/fdars-r/reference/sub-.irregFdata.md)
  : Subset method for irregFdata objects

## Random Processes

- [`r.bridge()`](https://sipemu.github.io/fdars-r/reference/r.bridge.md)
  : Generate Brownian Bridge
- [`r.brownian()`](https://sipemu.github.io/fdars-r/reference/r.brownian.md)
  : Generate Brownian Motion
- [`r.ou()`](https://sipemu.github.io/fdars-r/reference/r.ou.md) :
  Generate Ornstein-Uhlenbeck Process

## Statistical Tests

- [`fmean.test.fdata()`](https://sipemu.github.io/fdars-r/reference/fmean.test.fdata.md)
  : Test for Equality of Functional Means
- [`fequiv.test()`](https://sipemu.github.io/fdars-r/reference/fequiv.test.md)
  : Functional Equivalence Test (TOST)
- [`group.test()`](https://sipemu.github.io/fdars-r/reference/group.test.md)
  : Permutation Test for Group Differences

## Plotting

- [`autoplot(`*`<fdata>`*`)`](https://sipemu.github.io/fdars-r/reference/autoplot.fdata.md)
  : Create a ggplot for fdata objects
- [`plot(`*`<fdata>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.fdata.md)
  : Plot method for fdata objects
- [`boxplot(`*`<fdata>`*`)`](https://sipemu.github.io/fdars-r/reference/boxplot.fdata.md)
  : Functional Boxplot
- [`plot(`*`<fdata2pc>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.fdata2pc.md)
  : Plot FPCA Results
- [`plot(`*`<basis.auto>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.basis.auto.md)
  : Plot method for basis.auto objects
- [`plot(`*`<basis.cv>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.basis.cv.md)
  : Plot method for basis.cv objects
- [`plot(`*`<cluster.fcm>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.cluster.fcm.md)
  : Plot Method for cluster.fcm Objects
- [`plot(`*`<cluster.gmm>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.cluster.gmm.md)
  : Plot Method for cluster.gmm Objects
- [`plot(`*`<cluster.kmeans>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.cluster.kmeans.md)
  : Plot Method for cluster.kmeans Objects
- [`plot(`*`<cluster.optim>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.cluster.optim.md)
  : Plot Method for cluster.optim Objects
- [`plot(`*`<fanova>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.fanova.md)
  : Plot method for fanova objects
- [`plot(`*`<fclassif>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.fclassif.md)
  : Plot method for fclassif objects
- [`plot(`*`<fequiv.test>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.fequiv.test.md)
  : Plot method for fequiv.test
- [`plot(`*`<fmm>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.fmm.md)
  : Plot method for fmm objects
- [`plot(`*`<fosr>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.fosr.md)
  : Plot method for fosr objects
- [`plot(`*`<group.distance>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.group.distance.md)
  : Plot method for group.distance
- [`plot(`*`<outliergram>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.outliergram.md)
  : Plot Method for Outliergram Objects
- [`plot(`*`<outliers.fdata>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.outliers.fdata.md)
  : Plot method for outliers.fdata objects
- [`plot(`*`<pspline>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.pspline.md)
  : Plot method for pspline objects
- [`plot(`*`<pspline.2d>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.pspline.2d.md)
  : Plot method for pspline.2d objects
- [`plot(`*`<register.fd>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.register.fd.md)
  : Plot Method for register.fd Objects
- [`plot(`*`<magnitudeshape>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.magnitudeshape.md)
  : Plot Method for magnitudeshape Objects
- [`plot(`*`<amplitude_modulation>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.amplitude_modulation.md)
  : Plot method for amplitude_modulation objects
- [`plot(`*`<lomb_scargle_result>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.lomb_scargle_result.md)
  : Plot method for lomb_scargle_result objects
- [`plot(`*`<matrix_profile_result>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.matrix_profile_result.md)
  : Plot method for matrix_profile_result objects
- [`plot(`*`<peak_detection>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.peak_detection.md)
  : Plot method for peak_detection objects
- [`plot(`*`<peak_timing>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.peak_timing.md)
  : Plot method for peak_timing objects
- [`plot(`*`<ssa_result>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.ssa_result.md)
  : Plot method for ssa_result objects
- [`plot(`*`<stl_result>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.stl_result.md)
  : Plot method for stl_result objects

## Prediction

- [`predict(`*`<cluster.gmm>`*`)`](https://sipemu.github.io/fdars-r/reference/predict.cluster.gmm.md)
  : Predict Cluster Membership for New Functional Data
- [`predict(`*`<fosr>`*`)`](https://sipemu.github.io/fdars-r/reference/predict.fosr.md)
  : Predict from Function-on-Scalar Regression
- [`predict(`*`<fregre.fd>`*`)`](https://sipemu.github.io/fdars-r/reference/predict.fregre.fd.md)
  : Predict Method for Functional Regression (fregre.fd)
- [`predict(`*`<fregre.lm>`*`)`](https://sipemu.github.io/fdars-r/reference/predict.fregre.lm.md)
  : Predict method for fregre.lm objects
- [`predict(`*`<fregre.np>`*`)`](https://sipemu.github.io/fdars-r/reference/predict.fregre.np.md)
  : Predict Method for Nonparametric Functional Regression (fregre.np)
- [`predict(`*`<fregre.np.multi>`*`)`](https://sipemu.github.io/fdars-r/reference/predict.fregre.np.multi.md)
  : Predict method for fregre.np.multi

## Print & Summary

- [`print(`*`<fdata>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fdata.md)
  : Print method for fdata objects
- [`print(`*`<fdata2pc>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fdata2pc.md)
  : Print Method for FPCA Results
- [`print(`*`<fdata.bootstrap.ci>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fdata.bootstrap.ci.md)
  : Print method for bootstrap CI
- [`print(`*`<basis.auto>`*`)`](https://sipemu.github.io/fdars-r/reference/print.basis.auto.md)
  : Print method for basis.auto objects
- [`print(`*`<basis.cv>`*`)`](https://sipemu.github.io/fdars-r/reference/print.basis.cv.md)
  : Print method for basis.cv objects
- [`print(`*`<cluster.fcm>`*`)`](https://sipemu.github.io/fdars-r/reference/print.cluster.fcm.md)
  : Print Method for cluster.fcm Objects
- [`print(`*`<cluster.gmm>`*`)`](https://sipemu.github.io/fdars-r/reference/print.cluster.gmm.md)
  : Print Method for cluster.gmm Objects
- [`print(`*`<cluster.kmeans>`*`)`](https://sipemu.github.io/fdars-r/reference/print.cluster.kmeans.md)
  : Print Method for cluster.kmeans Objects
- [`print(`*`<cluster.optim>`*`)`](https://sipemu.github.io/fdars-r/reference/print.cluster.optim.md)
  : Print Method for cluster.optim Objects
- [`print(`*`<fanova>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fanova.md)
  : Print method for fanova objects
- [`print(`*`<fbplot>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fbplot.md)
  : Print Method for fbplot Objects
- [`print(`*`<fclassif>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fclassif.md)
  : Print method for fclassif objects
- [`print(`*`<fclassif.cv>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fclassif.cv.md)
  : Print method for fclassif.cv objects
- [`print(`*`<fmm>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fmm.md)
  : Print method for fmm objects
- [`print(`*`<fmm.test>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fmm.test.md)
  : Print method for fmm.test objects
- [`print(`*`<fosr>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fosr.md)
  : Print method for fosr objects
- [`print(`*`<fregre.fd>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fregre.fd.md)
  : Print method for fregre objects
- [`print(`*`<fregre.lm>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fregre.lm.md)
  : Print method for fregre.lm objects
- [`print(`*`<fregre.logistic>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fregre.logistic.md)
  : Print method for fregre.logistic objects
- [`print(`*`<fregre.np>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fregre.np.md)
  : Print method for fregre.np objects
- [`print(`*`<fregre.np.multi>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fregre.np.multi.md)
  : Print method for fregre.np.multi
- [`print(`*`<group.distance>`*`)`](https://sipemu.github.io/fdars-r/reference/print.group.distance.md)
  : Print method for group.distance
- [`print(`*`<fequiv.test>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fequiv.test.md)
  : Print method for fequiv.test
- [`print(`*`<group.test>`*`)`](https://sipemu.github.io/fdars-r/reference/print.group.test.md)
  : Print method for group.test
- [`print(`*`<kernel>`*`)`](https://sipemu.github.io/fdars-r/reference/print.kernel.md)
  : Print Method for Covariance Functions
- [`print(`*`<magnitudeshape>`*`)`](https://sipemu.github.io/fdars-r/reference/print.magnitudeshape.md)
  : Print Method for magnitudeshape Objects
- [`print(`*`<outliergram>`*`)`](https://sipemu.github.io/fdars-r/reference/print.outliergram.md)
  : Print Method for Outliergram Objects
- [`print(`*`<outliers.fdata>`*`)`](https://sipemu.github.io/fdars-r/reference/print.outliers.fdata.md)
  : Print method for outliers.fdata objects
- [`print(`*`<pspline>`*`)`](https://sipemu.github.io/fdars-r/reference/print.pspline.md)
  : Print method for pspline objects
- [`print(`*`<pspline.2d>`*`)`](https://sipemu.github.io/fdars-r/reference/print.pspline.2d.md)
  : Print method for pspline.2d objects
- [`print(`*`<register.fd>`*`)`](https://sipemu.github.io/fdars-r/reference/print.register.fd.md)
  : Print Method for register.fd Objects
- [`summary(`*`<basis.auto>`*`)`](https://sipemu.github.io/fdars-r/reference/summary.basis.auto.md)
  : Summary method for basis.auto objects
- [`summary(`*`<fdata>`*`)`](https://sipemu.github.io/fdars-r/reference/summary.fdata.md)
  : Summary method for fdata objects
- [`print(`*`<amplitude_modulation>`*`)`](https://sipemu.github.io/fdars-r/reference/print.amplitude_modulation.md)
  : Print method for amplitude_modulation objects
- [`print(`*`<autoperiod_result>`*`)`](https://sipemu.github.io/fdars-r/reference/print.autoperiod_result.md)
  : Print method for autoperiod_result objects
- [`print(`*`<cfd_autoperiod_result>`*`)`](https://sipemu.github.io/fdars-r/reference/print.cfd_autoperiod_result.md)
  : Print method for cfd_autoperiod_result objects
- [`print(`*`<decomposition>`*`)`](https://sipemu.github.io/fdars-r/reference/print.decomposition.md)
  : Print method for decomposition objects
- [`print(`*`<lomb_scargle_result>`*`)`](https://sipemu.github.io/fdars-r/reference/print.lomb_scargle_result.md)
  : Print method for lomb_scargle_result objects
- [`print(`*`<matrix_profile_result>`*`)`](https://sipemu.github.io/fdars-r/reference/print.matrix_profile_result.md)
  : Print method for matrix_profile_result objects
- [`print(`*`<multiFunData>`*`)`](https://sipemu.github.io/fdars-r/reference/print.multiFunData.md)
  : Print method for multiFunData objects
- [`print(`*`<multiple_periods>`*`)`](https://sipemu.github.io/fdars-r/reference/print.multiple_periods.md)
  : Print method for multiple_periods objects
- [`print(`*`<peak_detection>`*`)`](https://sipemu.github.io/fdars-r/reference/print.peak_detection.md)
  : Print method for peak_detection objects
- [`print(`*`<peak_timing>`*`)`](https://sipemu.github.io/fdars-r/reference/print.peak_timing.md)
  : Print method for peak_timing objects
- [`print(`*`<period_estimate>`*`)`](https://sipemu.github.io/fdars-r/reference/print.period_estimate.md)
  : Print method for period_estimate objects
- [`print(`*`<sazed_result>`*`)`](https://sipemu.github.io/fdars-r/reference/print.sazed_result.md)
  : Print method for sazed_result objects
- [`print(`*`<seasonality_changes>`*`)`](https://sipemu.github.io/fdars-r/reference/print.seasonality_changes.md)
  : Print method for seasonality_changes objects
- [`print(`*`<seasonality_changes_auto>`*`)`](https://sipemu.github.io/fdars-r/reference/print.seasonality_changes_auto.md)
  : Print method for seasonality_changes_auto objects
- [`print(`*`<seasonality_classification>`*`)`](https://sipemu.github.io/fdars-r/reference/print.seasonality_classification.md)
  : Print method for seasonality_classification objects
- [`print(`*`<ssa_result>`*`)`](https://sipemu.github.io/fdars-r/reference/print.ssa_result.md)
  : Print method for ssa_result objects
- [`print(`*`<stl_result>`*`)`](https://sipemu.github.io/fdars-r/reference/print.stl_result.md)
  : Print method for stl_result objects

## Other

- [`` `[`( ``*`<fdata>`*`)`](https://sipemu.github.io/fdars-r/reference/sub-.fdata.md)
  : Subset method for fdata objects
- [`Ops(`*`<fdata>`*`)`](https://sipemu.github.io/fdars-r/reference/Ops.fdata.md)
  : Arithmetic Operations for Functional Data
- [`kernels`](https://sipemu.github.io/fdars-r/reference/kernels.md) :
  Covariance Kernel Functions for Gaussian Processes
- [`fdars`](https://sipemu.github.io/fdars-r/reference/fdars-package.md)
  [`fdars-package`](https://sipemu.github.io/fdars-r/reference/fdars-package.md)
  : fdars: Functional Data Analysis in 'Rust'

## Internal

Low-level Rust wrappers and internal functions

- [`explain_ale_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_ale_logistic_rust.md)
  : ALE for logistic model
- [`explain_ale_rust()`](https://sipemu.github.io/fdars-r/reference/explain_ale_rust.md)
  : Accumulated Local Effects
- [`explain_anchor_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_anchor_logistic_rust.md)
  : Anchor for logistic model
- [`explain_anchor_rust()`](https://sipemu.github.io/fdars-r/reference/explain_anchor_rust.md)
  : Anchor explanation
- [`explain_beta_decomposition_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_beta_decomposition_logistic_rust.md)
  : Beta decomposition for logistic model
- [`explain_beta_decomposition_rust()`](https://sipemu.github.io/fdars-r/reference/explain_beta_decomposition_rust.md)
  : Beta decomposition by FPC components
- [`explain_calibration_rust()`](https://sipemu.github.io/fdars-r/reference/explain_calibration_rust.md)
  : Calibration diagnostics for logistic model
- [`explain_conditional_importance_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_conditional_importance_logistic_rust.md)
  : Conditional permutation importance for logistic model
- [`explain_conditional_importance_rust()`](https://sipemu.github.io/fdars-r/reference/explain_conditional_importance_rust.md)
  : Conditional permutation importance
- [`explain_conformal_rust()`](https://sipemu.github.io/fdars-r/reference/explain_conformal_rust.md)
  : Conformal prediction intervals
- [`explain_counterfactual_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_counterfactual_logistic_rust.md)
  : Counterfactual for logistic model (max_iter/step_size instead of
  target_value)
- [`explain_counterfactual_rust()`](https://sipemu.github.io/fdars-r/reference/explain_counterfactual_rust.md)
  : Counterfactual explanation
- [`explain_dfbetas_rust()`](https://sipemu.github.io/fdars-r/reference/explain_dfbetas_rust.md)
  : DFBETAS and DFFITS diagnostics
- [`explain_domain_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_domain_logistic_rust.md)
  : Domain selection for logistic model
- [`explain_domain_rust()`](https://sipemu.github.io/fdars-r/reference/explain_domain_rust.md)
  : Domain selection (important intervals)
- [`explain_ece_rust()`](https://sipemu.github.io/fdars-r/reference/explain_ece_rust.md)
  : Expected Calibration Error for logistic model
- [`explain_friedman_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_friedman_logistic_rust.md)
  : Friedman H-statistic for logistic model
- [`explain_friedman_rust()`](https://sipemu.github.io/fdars-r/reference/explain_friedman_rust.md)
  : Friedman H-statistic for interaction detection
- [`explain_importance_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_importance_logistic_rust.md)
  : Permutation importance for logistic model
- [`explain_importance_rust()`](https://sipemu.github.io/fdars-r/reference/explain_importance_rust.md)
  : FPC permutation importance
- [`explain_influence_rust()`](https://sipemu.github.io/fdars-r/reference/explain_influence_rust.md)
  : Influence diagnostics (leverage, Cook's distance)
- [`explain_lime_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_lime_logistic_rust.md)
  : LIME for logistic model
- [`explain_lime_rust()`](https://sipemu.github.io/fdars-r/reference/explain_lime_rust.md)
  : LIME local explanation
- [`explain_loo_rust()`](https://sipemu.github.io/fdars-r/reference/explain_loo_rust.md)
  : Leave-one-out cross-validation and PRESS
- [`explain_pdp_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_pdp_logistic_rust.md)
  : PDP for logistic model
- [`explain_pdp_rust()`](https://sipemu.github.io/fdars-r/reference/explain_pdp_rust.md)
  : Partial dependence plot for a single FPC component
- [`explain_pointwise_importance_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_pointwise_importance_logistic_rust.md)
  : Pointwise importance for logistic model
- [`explain_pointwise_importance_rust()`](https://sipemu.github.io/fdars-r/reference/explain_pointwise_importance_rust.md)
  : Pointwise importance of beta(t)
- [`explain_prediction_intervals_rust()`](https://sipemu.github.io/fdars-r/reference/explain_prediction_intervals_rust.md)
  : Prediction intervals for new observations
- [`explain_prototype_rust()`](https://sipemu.github.io/fdars-r/reference/explain_prototype_rust.md)
  : Prototype and criticism selection
- [`explain_regression_depth_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_regression_depth_logistic_rust.md)
  : Regression depth for logistic model
- [`explain_regression_depth_rust()`](https://sipemu.github.io/fdars-r/reference/explain_regression_depth_rust.md)
  : Regression depth
- [`explain_saliency_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_saliency_logistic_rust.md)
  : Saliency map for logistic model (fit-only, no data)
- [`explain_saliency_rust()`](https://sipemu.github.io/fdars-r/reference/explain_saliency_rust.md)
  : Functional saliency map
- [`explain_shap_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_shap_logistic_rust.md)
  : SHAP values for logistic model
- [`explain_shap_rust()`](https://sipemu.github.io/fdars-r/reference/explain_shap_rust.md)
  : SHAP values for FPC scores
- [`explain_significant_regions_rust()`](https://sipemu.github.io/fdars-r/reference/explain_significant_regions_rust.md)
  : Significant regions from beta(t) and its standard errors
- [`explain_sobol_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_sobol_logistic_rust.md)
  : Sobol indices for logistic model (no y, uses n_samples/seed)
- [`explain_sobol_rust()`](https://sipemu.github.io/fdars-r/reference/explain_sobol_rust.md)
  : Sobol sensitivity indices
- [`explain_stability_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_stability_logistic_rust.md)
  : Stability for logistic model (no fit param)
- [`explain_stability_rust()`](https://sipemu.github.io/fdars-r/reference/explain_stability_rust.md)
  : Explanation stability via bootstrap
- [`explain_vif_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/explain_vif_logistic_rust.md)
  : VIF for logistic model
- [`explain_vif_rust()`](https://sipemu.github.io/fdars-r/reference/explain_vif_rust.md)
  : Variance inflation factors for FPC scores
- [`elastic_amp_changepoint_rust()`](https://sipemu.github.io/fdars-r/reference/elastic_amp_changepoint_rust.md)
  : Amplitude changepoint detection
- [`elastic_fpca_changepoint_rust()`](https://sipemu.github.io/fdars-r/reference/elastic_fpca_changepoint_rust.md)
  : FPCA-based changepoint detection
- [`elastic_horiz_fpca_rust()`](https://sipemu.github.io/fdars-r/reference/elastic_horiz_fpca_rust.md)
  : Horizontal FPCA on warping functions
- [`elastic_joint_fpca_rust()`](https://sipemu.github.io/fdars-r/reference/elastic_joint_fpca_rust.md)
  : Joint FPCA (amplitude + phase)
- [`elastic_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/elastic_logistic_rust.md)
  : Elastic logistic regression
- [`elastic_pcr_attribution_rust()`](https://sipemu.github.io/fdars-r/reference/elastic_pcr_attribution_rust.md)
  : Elastic PCR attribution (amplitude vs phase importance)
- [`elastic_pcr_rust()`](https://sipemu.github.io/fdars-r/reference/elastic_pcr_rust.md)
  : Elastic principal component regression
- [`elastic_ph_changepoint_rust()`](https://sipemu.github.io/fdars-r/reference/elastic_ph_changepoint_rust.md)
  : Phase changepoint detection
- [`elastic_regression_rust()`](https://sipemu.github.io/fdars-r/reference/elastic_regression_rust.md)
  : Elastic regression (function-on-scalar with alignment)
- [`elastic_vert_fpca_rust()`](https://sipemu.github.io/fdars-r/reference/elastic_vert_fpca_rust.md)
  : Vertical FPCA on elastic-aligned data
- [`smooth_basis_gcv_rust()`](https://sipemu.github.io/fdars-r/reference/smooth_basis_gcv_rust.md)
  : Smooth functional data with GCV-selected penalty
- [`smooth_basis_rust()`](https://sipemu.github.io/fdars-r/reference/smooth_basis_rust.md)
  : Smooth functional data using B-spline or Fourier basis with fixed
  penalty
- [`bootstrap_ci_fregre_lm_rust()`](https://sipemu.github.io/fdars-r/reference/bootstrap_ci_fregre_lm_rust.md)
  : Bootstrap confidence intervals for β(t) in functional linear model
- [`bootstrap_ci_functional_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/bootstrap_ci_functional_logistic_rust.md)
  : Bootstrap confidence intervals for β(t) in functional logistic model
- [`conformal_classif_rust()`](https://sipemu.github.io/fdars-r/reference/conformal_classif_rust.md)
  : Conformal classification
- [`conformal_elastic_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/conformal_elastic_logistic_rust.md)
  : Conformal prediction for elastic logistic regression
- [`conformal_elastic_pcr_rust()`](https://sipemu.github.io/fdars-r/reference/conformal_elastic_pcr_rust.md)
  : Conformal prediction for elastic PCR
- [`conformal_elastic_regression_rust()`](https://sipemu.github.io/fdars-r/reference/conformal_elastic_regression_rust.md)
  : Conformal prediction for elastic regression
- [`conformal_fregre_lm_rust()`](https://sipemu.github.io/fdars-r/reference/conformal_fregre_lm_rust.md)
  : Conformal prediction for functional linear model
- [`conformal_fregre_np_rust()`](https://sipemu.github.io/fdars-r/reference/conformal_fregre_np_rust.md)
  : Conformal prediction for nonparametric functional regression
- [`conformal_generic_classification_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/conformal_generic_classification_logistic_rust.md)
  : Generic conformal classification for a pre-fitted logistic model
- [`conformal_generic_regression_lm_rust()`](https://sipemu.github.io/fdars-r/reference/conformal_generic_regression_lm_rust.md)
  : Generic conformal regression for a pre-fitted fregre.lm model
- [`conformal_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/conformal_logistic_rust.md)
  : Conformal prediction for logistic regression
- [`cv_conformal_classif_rust()`](https://sipemu.github.io/fdars-r/reference/cv_conformal_classif_rust.md)
  : CV-conformal classification
- [`cv_conformal_fregre_lm_rust()`](https://sipemu.github.io/fdars-r/reference/cv_conformal_fregre_lm_rust.md)
  : CV-conformal regression using functional linear model
- [`cv_conformal_fregre_np_rust()`](https://sipemu.github.io/fdars-r/reference/cv_conformal_fregre_np_rust.md)
  : CV-conformal regression using nonparametric kernel regression
- [`depth_rpd_1d_rust()`](https://sipemu.github.io/fdars-r/reference/depth_rpd_1d_rust.md)
  : Random Projection Depth with Derivatives (seeded)
- [`fanova_rust()`](https://sipemu.github.io/fdars-r/reference/fanova_rust.md)
  : Functional ANOVA
- [`fclassif_cv_rust()`](https://sipemu.github.io/fdars-r/reference/fclassif_cv_rust.md)
  : Cross-validated classification
- [`fclassif_dd_rust()`](https://sipemu.github.io/fdars-r/reference/fclassif_dd_rust.md)
  : DD-plot classification
- [`fclassif_kernel_rust()`](https://sipemu.github.io/fdars-r/reference/fclassif_kernel_rust.md)
  : Kernel classification
- [`fclassif_knn_rust()`](https://sipemu.github.io/fdars-r/reference/fclassif_knn_rust.md)
  : kNN classification
- [`fclassif_lda_rust()`](https://sipemu.github.io/fdars-r/reference/fclassif_lda_rust.md)
  : LDA classification
- [`fclassif_qda_rust()`](https://sipemu.github.io/fdars-r/reference/fclassif_qda_rust.md)
  : QDA classification
- [`fdata_gradient_rust()`](https://sipemu.github.io/fdars-r/reference/fdata_gradient_rust.md)
  : High-accuracy gradient using 5-point stencil (uniform) or 3-point
  Lagrange (non-uniform)
- [`fmm_predict_rust()`](https://sipemu.github.io/fdars-r/reference/fmm_predict_rust.md)
  : Predict from functional mixed model
- [`fmm_rust()`](https://sipemu.github.io/fdars-r/reference/fmm_rust.md)
  : Functional mixed model
- [`fmm_test_fixed_rust()`](https://sipemu.github.io/fdars-r/reference/fmm_test_fixed_rust.md)
  : Permutation test for fixed effects in FMM
- [`fosr_2d_rust()`](https://sipemu.github.io/fdars-r/reference/fosr_2d_rust.md)
  : 2D Function-on-scalar regression with tensor-product penalties
- [`fosr_fpc_rust()`](https://sipemu.github.io/fdars-r/reference/fosr_fpc_rust.md)
  : FPC-based function-on-scalar regression
- [`fosr_rust()`](https://sipemu.github.io/fdars-r/reference/fosr_rust.md)
  : Function-on-scalar regression (penalized)
- [`frcc_monitor_rust()`](https://sipemu.github.io/fdars-r/reference/frcc_monitor_rust.md)
  : FRCC Phase II: Monitor new data against a functional regression
  control chart.
- [`frcc_phase1_rust()`](https://sipemu.github.io/fdars-r/reference/frcc_phase1_rust.md)
  : FRCC Phase I: Build a functional regression control chart.
- [`fregre_cv_rust()`](https://sipemu.github.io/fdars-r/reference/fregre_cv_rust.md)
  : Cross-validation for FPC component selection
- [`fregre_huber_rust()`](https://sipemu.github.io/fdars-r/reference/fregre_huber_rust.md)
  : Huber M-estimation functional regression
- [`fregre_l1_rust()`](https://sipemu.github.io/fdars-r/reference/fregre_l1_rust.md)
  : L1 (LAD) functional regression
- [`fregre_lm_rust()`](https://sipemu.github.io/fdars-r/reference/fregre_lm_rust.md)
  : Functional linear model (FPC-based)
- [`fregre_np_mixed_rust()`](https://sipemu.github.io/fdars-r/reference/fregre_np_mixed_rust.md)
  : Nonparametric functional regression with mixed predictors
- [`functional_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/functional_logistic_rust.md)
  : Functional logistic regression
- [`gmm_cluster_rust()`](https://sipemu.github.io/fdars-r/reference/gmm_cluster_rust.md)
  : GMM clustering with automatic K selection
- [`gmm_em_rust()`](https://sipemu.github.io/fdars-r/reference/gmm_em_rust.md)
  : Raw GMM EM on feature matrix
- [`jackknife_plus_fregre_lm_rust()`](https://sipemu.github.io/fdars-r/reference/jackknife_plus_fregre_lm_rust.md)
  : Jackknife+ regression using functional linear model
- [`jackknife_plus_fregre_np_rust()`](https://sipemu.github.io/fdars-r/reference/jackknife_plus_fregre_np_rust.md)
  : Jackknife+ regression using nonparametric kernel regression
- [`mfpca_rust()`](https://sipemu.github.io/fdars-r/reference/mfpca_rust.md)
  : Multivariate FPCA.
- [`model_selection_ncomp_rust()`](https://sipemu.github.io/fdars-r/reference/model_selection_ncomp_rust.md)
  : Model selection for ncomp via AIC/BIC/GCV
- [`outliers_thres_lrt_with_dist_rust()`](https://sipemu.github.io/fdars-r/reference/outliers_thres_lrt_with_dist_rust.md)
  : LRT outlier threshold with full bootstrap distribution
- [`predict_elastic_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/predict_elastic_logistic_rust.md)
  : Predict from elastic logistic regression
- [`predict_elastic_regression_rust()`](https://sipemu.github.io/fdars-r/reference/predict_elastic_regression_rust.md)
  : Predict from elastic regression
- [`predict_fosr_2d_rust()`](https://sipemu.github.io/fdars-r/reference/predict_fosr_2d_rust.md)
  : Predict from 2D function-on-scalar regression
- [`predict_fosr_rust()`](https://sipemu.github.io/fdars-r/reference/predict_fosr_rust.md)
  : Predict from function-on-scalar regression
- [`predict_fregre_robust_rust()`](https://sipemu.github.io/fdars-r/reference/predict_fregre_robust_rust.md)
  : Predict from robust (L1/Huber) regression
- [`predict_functional_logistic_rust()`](https://sipemu.github.io/fdars-r/reference/predict_functional_logistic_rust.md)
  : Predict from functional logistic regression
- [`predict_gmm_rust()`](https://sipemu.github.io/fdars-r/reference/predict_gmm_rust.md)
  : Predict from GMM
- [`predict_scalar_on_shape_rust()`](https://sipemu.github.io/fdars-r/reference/predict_scalar_on_shape_rust.md)
  : Predict from a scalar-on-shape model
- [`scalar_on_shape_rust()`](https://sipemu.github.io/fdars-r/reference/scalar_on_shape_rust.md)
  : Fit scalar-on-shape regression model
- [`spm_ewma_rust()`](https://sipemu.github.io/fdars-r/reference/spm_ewma_rust.md)
  : EWMA-based SPM monitoring.
- [`spm_monitor_rust()`](https://sipemu.github.io/fdars-r/reference/spm_monitor_rust.md)
  : SPM Phase II: Monitor new data against an established chart.
- [`spm_phase1_rust()`](https://sipemu.github.io/fdars-r/reference/spm_phase1_rust.md)
  : SPM Phase I: Build a univariate control chart from in-control
  functional data.
- [`spm_spe_contrib_rust()`](https://sipemu.github.io/fdars-r/reference/spm_spe_contrib_rust.md)
  : SPE contribution diagnostics.
- [`spm_t2_contrib_rust()`](https://sipemu.github.io/fdars-r/reference/spm_t2_contrib_rust.md)
  : T-squared contribution diagnostics.
- [`tolerance_elastic_config_rust()`](https://sipemu.github.io/fdars-r/reference/tolerance_elastic_config_rust.md)
  : Elastic tolerance band with full config (amplitude + phase)
- [`tolerance_phase_rust()`](https://sipemu.github.io/fdars-r/reference/tolerance_phase_rust.md)
  : Phase tolerance band (warping variation)
- [`add_error_curve_1d()`](https://sipemu.github.io/fdars-r/reference/add_error_curve_1d.md)
  : Add curve-level Gaussian noise to functional data
- [`add_error_pointwise_1d()`](https://sipemu.github.io/fdars-r/reference/add_error_pointwise_1d.md)
  : Add pointwise Gaussian noise to functional data
- [`alignment_align_to_target()`](https://sipemu.github.io/fdars-r/reference/alignment_align_to_target.md)
  : Align all curves to a target curve
- [`alignment_amplitude_dist()`](https://sipemu.github.io/fdars-r/reference/alignment_amplitude_dist.md)
  : Amplitude self-distance matrix
- [`alignment_compose_warps()`](https://sipemu.github.io/fdars-r/reference/alignment_compose_warps.md)
  : Compose two warping functions
- [`alignment_constrained()`](https://sipemu.github.io/fdars-r/reference/alignment_constrained.md)
  : Elastic alignment with explicit landmark constraints
- [`alignment_cross_dist()`](https://sipemu.github.io/fdars-r/reference/alignment_cross_dist.md)
  : Elastic cross-distance matrix
- [`alignment_decomposition()`](https://sipemu.github.io/fdars-r/reference/alignment_decomposition.md)
  : Elastic phase-amplitude decomposition
- [`alignment_elastic_distance()`](https://sipemu.github.io/fdars-r/reference/alignment_elastic_distance.md)
  : Elastic (Fisher-Rao) distance between two curves
- [`alignment_elastic_pair()`](https://sipemu.github.io/fdars-r/reference/alignment_elastic_pair.md)
  : Elastic alignment of one curve to another
- [`alignment_karcher_mean()`](https://sipemu.github.io/fdars-r/reference/alignment_karcher_mean.md)
  : Karcher (Fréchet) mean in elastic metric
- [`alignment_pairwise_consistency()`](https://sipemu.github.io/fdars-r/reference/alignment_pairwise_consistency.md)
  : Pairwise alignment consistency
- [`alignment_phase_dist()`](https://sipemu.github.io/fdars-r/reference/alignment_phase_dist.md)
  : Phase self-distance matrix
- [`alignment_quality_compute()`](https://sipemu.github.io/fdars-r/reference/alignment_quality_compute.md)
  : Compute alignment quality metrics
- [`alignment_reparameterize()`](https://sipemu.github.io/fdars-r/reference/alignment_reparameterize.md)
  : Apply warping function to reparameterize a curve
- [`alignment_self_dist()`](https://sipemu.github.io/fdars-r/reference/alignment_self_dist.md)
  : Elastic self-distance matrix
- [`alignment_srsf_inverse()`](https://sipemu.github.io/fdars-r/reference/alignment_srsf_inverse.md)
  : Inverse SRSF: reconstruct curve from SRSF representation
- [`alignment_srsf_transform()`](https://sipemu.github.io/fdars-r/reference/alignment_srsf_transform.md)
  : SRSF transform of functional data
- [`alignment_tsrvf_from_karcher()`](https://sipemu.github.io/fdars-r/reference/alignment_tsrvf_from_karcher.md)
  : Compute TSRVF from a pre-computed Karcher mean
- [`alignment_tsrvf_inverse()`](https://sipemu.github.io/fdars-r/reference/alignment_tsrvf_inverse.md)
  : Inverse TSRVF: reconstruct curves from tangent vectors
- [`alignment_tsrvf_transform()`](https://sipemu.github.io/fdars-r/reference/alignment_tsrvf_transform.md)
  : Full TSRVF transform: compute Karcher mean + transport to tangent
  space
- [`alignment_warp_complexity()`](https://sipemu.github.io/fdars-r/reference/alignment_warp_complexity.md)
  : Compute warp complexity (geodesic distance from identity)
- [`alignment_warp_smoothness()`](https://sipemu.github.io/fdars-r/reference/alignment_warp_smoothness.md)
  : Compute warp smoothness (bending energy)
- [`alignment_with_landmarks()`](https://sipemu.github.io/fdars-r/reference/alignment_with_landmarks.md)
  : Elastic alignment with automatic landmark detection
- [`andrews_loadings()`](https://sipemu.github.io/fdars-r/reference/andrews_loadings.md)
  : Andrews Loadings: Project FPCA Eigenfunctions to Original Variables
- [`andrews_transform()`](https://sipemu.github.io/fdars-r/reference/andrews_transform.md)
  : Andrews Transformation
- [`basis_aic_1d()`](https://sipemu.github.io/fdars-r/reference/basis_aic_1d.md)
  : Compute AIC for basis fit AIC = n \* log(RSS/n) + 2 \* total_edf
  Where total_edf = n_curves \* edf (each curve has edf parameters) When
  pooled=true: compute single AIC across all curves When pooled=false:
  compute per-curve AIC and return mean
- [`basis_bic_1d()`](https://sipemu.github.io/fdars-r/reference/basis_bic_1d.md)
  : Compute BIC for basis fit BIC = n \* log(RSS/n) + log(n) \*
  total_edf Where total_edf = n_curves \* edf (each curve has edf
  parameters) When pooled=true: compute single BIC across all curves
  When pooled=false: compute per-curve BIC and return mean
- [`basis_gcv_1d()`](https://sipemu.github.io/fdars-r/reference/basis_gcv_1d.md)
  : Compute GCV score for basis fit GCV = RSS/n / (1 - edf/n)^2 When
  pooled=true: compute single GCV across all curves When pooled=false:
  compute per-curve GCV and return mean
- [`calinski_harabasz()`](https://sipemu.github.io/fdars-r/reference/calinski_harabasz.md)
  : Compute Calinski-Harabasz index (variance ratio criterion) Higher
  values indicate better defined clusters
- [`compute_adot()`](https://sipemu.github.io/fdars-r/reference/compute_adot.md)
  : Compute the Adot matrix (parallelized)
- [`depth_bd_1d()`](https://sipemu.github.io/fdars-r/reference/depth_bd_1d.md)
  : Band Depth (BD) for 1D functional data BD(x) = proportion of pairs
  (i,j) where x lies within the band formed by curves i and j A curve
  lies in the band if at every time point t, min(X_i(t), X_j(t)) \<=
  x(t) \<= max(X_i(t), X_j(t))
- [`depth_fm_1d()`](https://sipemu.github.io/fdars-r/reference/depth_fm_1d.md)
  : Compute Fraiman-Muniz depth
- [`depth_fm_2d()`](https://sipemu.github.io/fdars-r/reference/depth_fm_2d.md)
  : Fraiman-Muniz depth for 2D functional data (surfaces) Integrates
  univariate depth over (s,t) grid
- [`depth_fsd_1d()`](https://sipemu.github.io/fdars-r/reference/depth_fsd_1d.md)
  : Compute Functional Spatial Depth
- [`depth_fsd_2d()`](https://sipemu.github.io/fdars-r/reference/depth_fsd_2d.md)
  : Functional Spatial Depth for 2D functional data
- [`depth_kfsd_1d()`](https://sipemu.github.io/fdars-r/reference/depth_kfsd_1d.md)
  : Kernel Functional Spatial Depth (KFSD) for 1D functional data
  Implements the RKHS-based formulation matching fda.usc h is treated as
  the actual bandwidth, matching how fda.usc uses hq2 argvals is used
  for trapezoidal integration to compute L2 norms
- [`depth_kfsd_2d()`](https://sipemu.github.io/fdars-r/reference/depth_kfsd_2d.md)
  : Kernel Functional Spatial Depth (KFSD) for 2D functional data
  Implements the RKHS-based formulation matching fda.usc
- [`depth_mbd_1d()`](https://sipemu.github.io/fdars-r/reference/depth_mbd_1d.md)
  : Modified Band Depth (MBD) for 1D functional data MBD(x) = average
  over pairs (i,j) of the proportion of the domain where x is inside the
  band This is more robust than BD as it doesn't require complete
  containment
- [`depth_mei_1d()`](https://sipemu.github.io/fdars-r/reference/depth_mei_1d.md)
  : Modified Epigraph Index (MEI) for 1D functional data MEI measures
  the proportion of time a curve is below other curves MEI(x_i) = (1/n)
  \* sum_j (1/m) \* sum_t I(x_i(t) \< x_j(t)) + 0.5\*I(x_i(t) = x_j(t))
- [`depth_mode_1d()`](https://sipemu.github.io/fdars-r/reference/depth_mode_1d.md)
  : Compute modal depth
- [`depth_mode_2d()`](https://sipemu.github.io/fdars-r/reference/depth_mode_2d.md)
  : Modal depth for 2D functional data (surfaces) Uses L2 distance in
  the flattened surface space
- [`depth_rp_1d()`](https://sipemu.github.io/fdars-r/reference/depth_rp_1d.md)
  : Compute random projection depth
- [`depth_rp_1d_seeded()`](https://sipemu.github.io/fdars-r/reference/depth_rp_1d_seeded.md)
  : Random projection depth with optional seed
- [`depth_rp_2d()`](https://sipemu.github.io/fdars-r/reference/depth_rp_2d.md)
  : Random projection depth for 2D functional data (surfaces) Projects
  surfaces to scalars using random projections
- [`depth_rt_1d()`](https://sipemu.github.io/fdars-r/reference/depth_rt_1d.md)
  : Compute random Tukey depth
- [`depth_rt_1d_seeded()`](https://sipemu.github.io/fdars-r/reference/depth_rt_1d_seeded.md)
  : Random Tukey depth with optional seed
- [`depth_rt_2d()`](https://sipemu.github.io/fdars-r/reference/depth_rt_2d.md)
  : Random Tukey depth for 2D functional data (surfaces)
- [`detect_amplitude_modulation()`](https://sipemu.github.io/fdars-r/reference/detect_amplitude_modulation.md)
  : Detect Amplitude Modulation in Seasonal Time Series
- [`df_to_fdata2d()`](https://sipemu.github.io/fdars-r/reference/df_to_fdata2d.md)
  : Convert DataFrame to 2D functional data
- [`fdata2basis_cv()`](https://sipemu.github.io/fdars-r/reference/fdata2basis_cv.md)
  : Cross-Validation for Basis Function Number Selection
- [`fdata_center_1d()`](https://sipemu.github.io/fdars-r/reference/fdata_center_1d.md)
  : Center functional data by subtracting the mean function
- [`fdata_deriv_1d()`](https://sipemu.github.io/fdars-r/reference/fdata_deriv_1d.md)
  : Compute numerical derivative of functional data (parallelized over
  rows)
- [`fdata_deriv_2d()`](https://sipemu.github.io/fdars-r/reference/fdata_deriv_2d.md)
  : Compute 2D partial derivatives for surface data
- [`fdata_mean_1d()`](https://sipemu.github.io/fdars-r/reference/fdata_mean_1d.md)
  : Compute the mean function across all samples (1D)
- [`fdata_mean_2d()`](https://sipemu.github.io/fdars-r/reference/fdata_mean_2d.md)
  : Compute the mean function across all samples (2D surfaces) Data is
  stored as n x (m1\*m2) matrix where each row is a flattened surface
- [`fdata_norm_lp_1d()`](https://sipemu.github.io/fdars-r/reference/fdata_norm_lp_1d.md)
  : Compute Lp norm for each sample
- [`fuzzycmeans_fd()`](https://sipemu.github.io/fdars-r/reference/fuzzycmeans_fd.md)
  : Fuzzy C-Means clustering for functional data m_fuzz is the fuzziness
  parameter (typically 2)
- [`geometric_median_1d()`](https://sipemu.github.io/fdars-r/reference/geometric_median_1d.md)
  : Compute the geometric median (L1 median) of functional data using
  Weiszfeld's algorithm The geometric median minimizes sum of L2
  distances to all curves
- [`geometric_median_2d()`](https://sipemu.github.io/fdars-r/reference/geometric_median_2d.md)
  : Compute the geometric median (L1 median) of 2D functional data using
  Weiszfeld's algorithm Data is stored as n x (m1\*m2) matrix where each
  row is a flattened surface
- [`inprod_fdata()`](https://sipemu.github.io/fdars-r/reference/inprod_fdata.md)
  : Inner product of two functional data objects \<f, g\> =
  integral(f(t) \* g(t) dt)
- [`int_simpson()`](https://sipemu.github.io/fdars-r/reference/int_simpson.md)
  : Simpson's rule integration for functional data Integrates each curve
  over the domain
- [`irreg_fdata2basis()`](https://sipemu.github.io/fdars-r/reference/irreg_fdata2basis.md)
  : Fit basis functions to irregular functional data Each curve is
  individually fitted via least squares at its own observation points
  basis_type: 0 = bspline, 1 = fourier
- [`irreg_integrate()`](https://sipemu.github.io/fdars-r/reference/irreg_integrate.md)
  : Compute integral for each curve in irregular functional data
- [`irreg_mean_kernel()`](https://sipemu.github.io/fdars-r/reference/irreg_mean_kernel.md)
  : Estimate mean function for irregular data using kernel smoothing
- [`irreg_metric_lp()`](https://sipemu.github.io/fdars-r/reference/irreg_metric_lp.md)
  : Compute pairwise Lp distances for irregular functional data
- [`irreg_norm_lp()`](https://sipemu.github.io/fdars-r/reference/irreg_norm_lp.md)
  : Compute Lp norm for each curve in irregular functional data
- [`irreg_to_regular()`](https://sipemu.github.io/fdars-r/reference/irreg_to_regular.md)
  : Convert irregular data to regular grid via interpolation
- [`kmeans_fd()`](https://sipemu.github.io/fdars-r/reference/kmeans_fd.md)
  : Functional k-means clustering
- [`knn_gcv()`](https://sipemu.github.io/fdars-r/reference/knn_gcv.md) :
  k-NN with Global Cross-Validation Finds a single optimal k for all
  observations
- [`knn_lcv()`](https://sipemu.github.io/fdars-r/reference/knn_lcv.md) :
  k-NN with Local Cross-Validation Finds an optimal k for each
  observation
- [`knn_predict()`](https://sipemu.github.io/fdars-r/reference/knn_predict.md)
  : Kernel prediction with fixed bandwidth for prediction on new data
- [`landmark_detect()`](https://sipemu.github.io/fdars-r/reference/landmark_detect.md)
  : Detect landmarks in a single curve
- [`landmark_register_curves()`](https://sipemu.github.io/fdars-r/reference/landmark_register_curves.md)
  : Detect landmarks and register curves
- [`metric_dtw_cross_1d()`](https://sipemu.github.io/fdars-r/reference/metric_dtw_cross_1d.md)
  : Compute DTW distance matrix for cross-distances (n1 x n2)
- [`metric_dtw_self_1d()`](https://sipemu.github.io/fdars-r/reference/metric_dtw_self_1d.md)
  : Compute DTW distance matrix for self-distances (symmetric)
- [`metric_hausdorff_1d()`](https://sipemu.github.io/fdars-r/reference/metric_hausdorff_1d.md)
  : Compute Hausdorff distance matrix for self-distances (symmetric)
- [`metric_hausdorff_2d()`](https://sipemu.github.io/fdars-r/reference/metric_hausdorff_2d.md)
  : Compute Hausdorff distance for 2D functional data (surfaces)
- [`metric_hausdorff_cross_1d()`](https://sipemu.github.io/fdars-r/reference/metric_hausdorff_cross_1d.md)
  : Compute Hausdorff distance matrix for cross-distances (n1 x n2)
- [`metric_hausdorff_cross_2d()`](https://sipemu.github.io/fdars-r/reference/metric_hausdorff_cross_2d.md)
  : Compute Hausdorff cross-distances for 2D functional data
- [`metric_kl_cross_1d()`](https://sipemu.github.io/fdars-r/reference/metric_kl_cross_1d.md)
  : Compute symmetric KL divergence matrix for cross-distances (1D)
- [`metric_kl_self_1d()`](https://sipemu.github.io/fdars-r/reference/metric_kl_self_1d.md)
  : Compute symmetric KL divergence matrix for self-distances (1D)
  Curves are first normalized to be valid probability distributions
- [`metric_lp_1d()`](https://sipemu.github.io/fdars-r/reference/metric_lp_1d.md)
  : Compute Lp distance matrix between two sets of functional data
- [`metric_lp_2d()`](https://sipemu.github.io/fdars-r/reference/metric_lp_2d.md)
  : Compute Lp distance between two 2D functional data objects
  (surfaces)
- [`metric_lp_self_1d()`](https://sipemu.github.io/fdars-r/reference/metric_lp_self_1d.md)
  : Compute Lp distance matrix for self-distances (symmetric)
- [`metric_lp_self_2d()`](https://sipemu.github.io/fdars-r/reference/metric_lp_self_2d.md)
  : Compute Lp self-distance matrix for 2D functional data (symmetric)
- [`metric_soft_dtw_barycenter()`](https://sipemu.github.io/fdars-r/reference/metric_soft_dtw_barycenter.md)
  : Soft-DTW barycenter computation
- [`metric_soft_dtw_cross_1d()`](https://sipemu.github.io/fdars-r/reference/metric_soft_dtw_cross_1d.md)
  : Soft-DTW cross-distance matrix
- [`metric_soft_dtw_div_cross_1d()`](https://sipemu.github.io/fdars-r/reference/metric_soft_dtw_div_cross_1d.md)
  : Soft-DTW divergence cross-distance matrix
- [`metric_soft_dtw_div_self_1d()`](https://sipemu.github.io/fdars-r/reference/metric_soft_dtw_div_self_1d.md)
  : Soft-DTW divergence self-distance matrix
- [`metric_soft_dtw_self_1d()`](https://sipemu.github.io/fdars-r/reference/metric_soft_dtw_self_1d.md)
  : Soft-DTW self-distance matrix
- [`outlier_summary()`](https://sipemu.github.io/fdars-r/reference/outlier_summary.md)
  : Unified Outlier Summary
- [`outliers_lrt()`](https://sipemu.github.io/fdars-r/reference/outliers_lrt.md)
  : LRT-based outlier detection Returns indices of detected outliers
- [`outliers_thres_lrt()`](https://sipemu.github.io/fdars-r/reference/outliers_thres_lrt.md)
  : Compute bootstrap threshold for LRT outlier detection Highly
  parallelized across bootstrap iterations
- [`pcvm_statistic()`](https://sipemu.github.io/fdars-r/reference/pcvm_statistic.md)
  : Compute the PCvM statistic
- [`plot(`*`<amplitude_modulation>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.amplitude_modulation.md)
  : Plot method for amplitude_modulation objects
- [`plot(`*`<lomb_scargle_result>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.lomb_scargle_result.md)
  : Plot method for lomb_scargle_result objects
- [`plot(`*`<matrix_profile_result>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.matrix_profile_result.md)
  : Plot method for matrix_profile_result objects
- [`plot(`*`<peak_detection>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.peak_detection.md)
  : Plot method for peak_detection objects
- [`plot(`*`<peak_timing>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.peak_timing.md)
  : Plot method for peak_timing objects
- [`plot(`*`<ssa_result>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.ssa_result.md)
  : Plot method for ssa_result objects
- [`plot(`*`<stl_result>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.stl_result.md)
  : Plot method for stl_result objects
- [`print(`*`<amplitude_modulation>`*`)`](https://sipemu.github.io/fdars-r/reference/print.amplitude_modulation.md)
  : Print method for amplitude_modulation objects
- [`print(`*`<autoperiod_result>`*`)`](https://sipemu.github.io/fdars-r/reference/print.autoperiod_result.md)
  : Print method for autoperiod_result objects
- [`print(`*`<cfd_autoperiod_result>`*`)`](https://sipemu.github.io/fdars-r/reference/print.cfd_autoperiod_result.md)
  : Print method for cfd_autoperiod_result objects
- [`print(`*`<lomb_scargle_result>`*`)`](https://sipemu.github.io/fdars-r/reference/print.lomb_scargle_result.md)
  : Print method for lomb_scargle_result objects
- [`print(`*`<matrix_profile_result>`*`)`](https://sipemu.github.io/fdars-r/reference/print.matrix_profile_result.md)
  : Print method for matrix_profile_result objects
- [`print(`*`<multiple_periods>`*`)`](https://sipemu.github.io/fdars-r/reference/print.multiple_periods.md)
  : Print method for multiple_periods objects
- [`print(`*`<peak_detection>`*`)`](https://sipemu.github.io/fdars-r/reference/print.peak_detection.md)
  : Print method for peak_detection objects
- [`print(`*`<peak_timing>`*`)`](https://sipemu.github.io/fdars-r/reference/print.peak_timing.md)
  : Print method for peak_timing objects
- [`print(`*`<period_estimate>`*`)`](https://sipemu.github.io/fdars-r/reference/print.period_estimate.md)
  : Print method for period_estimate objects
- [`print(`*`<sazed_result>`*`)`](https://sipemu.github.io/fdars-r/reference/print.sazed_result.md)
  : Print method for sazed_result objects
- [`print(`*`<seasonality_changes>`*`)`](https://sipemu.github.io/fdars-r/reference/print.seasonality_changes.md)
  : Print method for seasonality_changes objects
- [`print(`*`<seasonality_changes_auto>`*`)`](https://sipemu.github.io/fdars-r/reference/print.seasonality_changes_auto.md)
  : Print method for seasonality_changes_auto objects
- [`print(`*`<seasonality_classification>`*`)`](https://sipemu.github.io/fdars-r/reference/print.seasonality_classification.md)
  : Print method for seasonality_classification objects
- [`print(`*`<ssa_result>`*`)`](https://sipemu.github.io/fdars-r/reference/print.ssa_result.md)
  : Print method for ssa_result objects
- [`print(`*`<stl_result>`*`)`](https://sipemu.github.io/fdars-r/reference/print.stl_result.md)
  : Print method for stl_result objects
- [`pspline_fit_1d()`](https://sipemu.github.io/fdars-r/reference/pspline_fit_1d.md)
  : P-spline fitting: returns coefficients, fitted values, and
  diagnostics
- [`pspline_fit_2d()`](https://sipemu.github.io/fdars-r/reference/pspline_fit_2d.md)
  : 2D P-spline fitting with anisotropic penalties
- [`register_shift_1d()`](https://sipemu.github.io/fdars-r/reference/register_shift_1d.md)
  : Shift registration: find optimal horizontal shift for each curve to
  align with a target (usually the mean)
- [`rp_stat()`](https://sipemu.github.io/fdars-r/reference/rp_stat.md) :
  Compute random projection statistics (parallelized over projections)
- [`s_knn()`](https://sipemu.github.io/fdars-r/reference/s_knn.md) :
  K-Nearest Neighbors smoother matrix
- [`s_llr()`](https://sipemu.github.io/fdars-r/reference/s_llr.md) :
  Local Linear Regression smoother matrix Uses weighted least squares
  with degree-1 polynomial
- [`s_lpr()`](https://sipemu.github.io/fdars-r/reference/s_lpr.md) :
  Local Polynomial Regression smoother matrix Solves (p+1)×(p+1)
  weighted least squares system for each point
- [`s_nw()`](https://sipemu.github.io/fdars-r/reference/s_nw.md) :
  Nadaraya-Watson smoother matrix S_ij = K((t_i - t_j)/h) \* w_j /
  sum_k(K((t_i - t_k)/h) \* w_k)
- [`scale_minmax()`](https://sipemu.github.io/fdars-r/reference/scale_minmax.md)
  : Min-Max scaling for functional data
- [`seasonal_analyze_peak_timing()`](https://sipemu.github.io/fdars-r/reference/seasonal_analyze_peak_timing.md)
  : Analyze peak timing variability across cycles (uses Fourier
  smoothing)
- [`seasonal_autoperiod()`](https://sipemu.github.io/fdars-r/reference/seasonal_autoperiod.md)
  : Autoperiod: Hybrid FFT + ACF period detection with gradient ascent
  refinement Returns period, confidence, FFT power, ACF validation
  score, and candidates
- [`seasonal_cfd_autoperiod()`](https://sipemu.github.io/fdars-r/reference/seasonal_cfd_autoperiod.md)
  : CFDAutoperiod: Clustered Filtered Detrended Autoperiod Uses
  differencing for detrending and clustering for robust period detection
- [`seasonal_classify_seasonality()`](https://sipemu.github.io/fdars-r/reference/seasonal_classify_seasonality.md)
  : Classify seasonality type
- [`seasonal_decompose()`](https://sipemu.github.io/fdars-r/reference/seasonal_decompose.md)
  : Decompose functional data into trend, seasonal, and remainder
  components
- [`seasonal_detect_amplitude_modulation()`](https://sipemu.github.io/fdars-r/reference/seasonal_detect_amplitude_modulation.md)
  : Detect amplitude modulation in seasonal time series using Hilbert
  transform
- [`seasonal_detect_amplitude_modulation_wavelet()`](https://sipemu.github.io/fdars-r/reference/seasonal_detect_amplitude_modulation_wavelet.md)
  : Detect amplitude modulation using wavelet transform (Morlet wavelet)
- [`seasonal_detect_changes()`](https://sipemu.github.io/fdars-r/reference/seasonal_detect_changes.md)
  : Detect seasonality changes (onset/cessation)
- [`seasonal_detect_changes_auto()`](https://sipemu.github.io/fdars-r/reference/seasonal_detect_changes_auto.md)
  : Detect seasonality changes with automatic threshold
- [`seasonal_detect_multiple_periods()`](https://sipemu.github.io/fdars-r/reference/seasonal_detect_multiple_periods.md)
  : Detect multiple concurrent periodicities using iterative residual
  subtraction
- [`seasonal_detect_peaks()`](https://sipemu.github.io/fdars-r/reference/seasonal_detect_peaks.md)
  : Detect peaks in functional data using Fourier basis smoothing
- [`seasonal_detrend()`](https://sipemu.github.io/fdars-r/reference/seasonal_detrend.md)
  : Detrend functional data using specified method Returns trend,
  detrended data, method used, RSS per curve, and number of parameters
- [`seasonal_estimate_period_acf()`](https://sipemu.github.io/fdars-r/reference/seasonal_estimate_period_acf.md)
  : Estimate period using autocorrelation
- [`seasonal_estimate_period_fft()`](https://sipemu.github.io/fdars-r/reference/seasonal_estimate_period_fft.md)
  : Estimate period using FFT periodogram
- [`seasonal_instantaneous_period()`](https://sipemu.github.io/fdars-r/reference/seasonal_instantaneous_period.md)
  : Estimate instantaneous period using Hilbert transform
- [`seasonal_lomb_scargle()`](https://sipemu.github.io/fdars-r/reference/seasonal_lomb_scargle.md)
  : Lomb-Scargle periodogram for irregularly sampled data Computes the
  power spectrum and significance for period detection
- [`seasonal_matrix_profile()`](https://sipemu.github.io/fdars-r/reference/seasonal_matrix_profile.md)
  : Matrix Profile for motif discovery and period detection Uses STOMP
  algorithm for efficient computation
- [`seasonal_sazed()`](https://sipemu.github.io/fdars-r/reference/seasonal_sazed.md)
  : SAZED: Spectral-ACF Zero-crossing Ensemble Detection A
  parameter-free ensemble method for robust period detection Returns
  period, confidence, component periods, and agreeing component count
- [`seasonal_ssa()`](https://sipemu.github.io/fdars-r/reference/seasonal_ssa.md)
  : Singular Spectrum Analysis for time series decomposition Extracts
  trend, seasonal, and noise components via SVD
- [`seasonal_stl()`](https://sipemu.github.io/fdars-r/reference/seasonal_stl.md)
  : STL (Seasonal and Trend decomposition using LOESS) Implements
  Cleveland et al. 1990 algorithm
- [`seasonal_strength_spectral()`](https://sipemu.github.io/fdars-r/reference/seasonal_strength_spectral.md)
  : Measure seasonal strength using spectral method
- [`seasonal_strength_variance()`](https://sipemu.github.io/fdars-r/reference/seasonal_strength_variance.md)
  : Measure seasonal strength using variance decomposition
- [`seasonal_strength_wavelet()`](https://sipemu.github.io/fdars-r/reference/seasonal_strength_wavelet.md)
  : Measure seasonal strength using wavelet (Morlet) method
- [`seasonal_strength_windowed()`](https://sipemu.github.io/fdars-r/reference/seasonal_strength_windowed.md)
  : Time-varying seasonal strength using sliding windows
- [`select_basis_auto()`](https://sipemu.github.io/fdars-r/reference/select_basis_auto.md)
  : Automatic basis selection for each curve individually.
- [`semimetric_basis_cross_1d()`](https://sipemu.github.io/fdars-r/reference/semimetric_basis_cross_1d.md)
  : Basis coefficient semimetric (cross-distances)
- [`semimetric_basis_self_1d()`](https://sipemu.github.io/fdars-r/reference/semimetric_basis_self_1d.md)
  : Basis coefficient semimetric (self-distances)
- [`semimetric_deriv_cross_1d()`](https://sipemu.github.io/fdars-r/reference/semimetric_deriv_cross_1d.md)
  : Derivative-based semimetric (cross-distances)
- [`semimetric_deriv_self_1d()`](https://sipemu.github.io/fdars-r/reference/semimetric_deriv_self_1d.md)
  : Derivative-based semimetric (self-distances)
- [`semimetric_fourier_cross_1d()`](https://sipemu.github.io/fdars-r/reference/semimetric_fourier_cross_1d.md)
  : Compute semimetric based on Fourier coefficients for cross-distances
- [`semimetric_fourier_self_1d()`](https://sipemu.github.io/fdars-r/reference/semimetric_fourier_self_1d.md)
  : Compute semimetric based on Fourier coefficients for self-distances
  (symmetric) Uses FFT to compute Fourier coefficients and then L2
  distance on coefficients
- [`semimetric_hshift_cross_1d()`](https://sipemu.github.io/fdars-r/reference/semimetric_hshift_cross_1d.md)
  : Compute semimetric based on horizontal shift for cross-distances
- [`semimetric_hshift_self_1d()`](https://sipemu.github.io/fdars-r/reference/semimetric_hshift_self_1d.md)
  : Compute semimetric based on horizontal shift for self-distances
  (symmetric) This finds the minimum L2 distance after optimally
  shifting one curve horizontally
- [`semimetric_pca_cross_1d()`](https://sipemu.github.io/fdars-r/reference/semimetric_pca_cross_1d.md)
  : PCA-based semimetric (cross-distances)
- [`semimetric_pca_self_1d()`](https://sipemu.github.io/fdars-r/reference/semimetric_pca_self_1d.md)
  : PCA-based semimetric (self-distances)
- [`silhouette_score()`](https://sipemu.github.io/fdars-r/reference/silhouette_score.md)
  : Compute silhouette score for clustering Returns the mean silhouette
  coefficient across all samples
- [`sim_kl_1d()`](https://sipemu.github.io/fdars-r/reference/sim_kl_1d.md)
  : Simulate functional data via Karhunen-Loève expansion
- [`streaming_depth_batch()`](https://sipemu.github.io/fdars-r/reference/streaming_depth_batch.md)
  : Streaming depth: batch self-depth computation
- [`streaming_depth_one()`](https://sipemu.github.io/fdars-r/reference/streaming_depth_one.md)
  : Streaming depth: single curve against reference
- [`streaming_depth_vs_ref()`](https://sipemu.github.io/fdars-r/reference/streaming_depth_vs_ref.md)
  : Streaming depth: new data against reference
- [`tolerance_conformal()`](https://sipemu.github.io/fdars-r/reference/tolerance_conformal.md)
  : Conformal prediction band
- [`tolerance_elastic()`](https://sipemu.github.io/fdars-r/reference/tolerance_elastic.md)
  : Elastic tolerance band (alignment + FPCA)
- [`tolerance_exponential()`](https://sipemu.github.io/fdars-r/reference/tolerance_exponential.md)
  : Exponential family tolerance band
- [`tolerance_fpca()`](https://sipemu.github.io/fdars-r/reference/tolerance_fpca.md)
  : FPCA-based tolerance band
- [`tolerance_scb_degras()`](https://sipemu.github.io/fdars-r/reference/tolerance_scb_degras.md)
  : SCB mean confidence band (Degras method)
- [`basis2fdata_1d()`](https://sipemu.github.io/fdars-r/reference/basis2fdata_1d.md)
  : Reconstruct functional data from basis coefficients Returns data
  matrix (n x m)
- [`eigenfunctions_1d()`](https://sipemu.github.io/fdars-r/reference/eigenfunctions_1d.md)
  : Compute eigenfunction basis values efun_type: 0 = Fourier, 1 = Poly,
  2 = PolyHigh, 3 = Wiener
- [`eigenvalues_1d()`](https://sipemu.github.io/fdars-r/reference/eigenvalues_1d.md)
  : Generate eigenvalue sequence eval_type: 0 = linear, 1 = exponential,
  2 = wiener
- [`fdata2basis_1d()`](https://sipemu.github.io/fdars-r/reference/fdata2basis_1d.md)
  : Convert functional data to basis coefficients type: 0 = bspline, 1 =
  fourier
- [`fdata2pc_1d()`](https://sipemu.github.io/fdars-r/reference/fdata2pc_1d.md)
  : Perform functional PCA via SVD on centered data Returns: singular
  values, rotation matrix (loadings), scores, mean
- [`fdata2pls_1d()`](https://sipemu.github.io/fdars-r/reference/fdata2pls_1d.md)
  : Perform PLS via NIPALS algorithm Returns: weights, scores, loadings
