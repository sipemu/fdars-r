# Package index

## Functional Data Objects

- [`fdata()`](https://sipemu.github.io/fdars-r/reference/fdata.md) :
  Create a functional data object
- [`fdata.cen()`](https://sipemu.github.io/fdars-r/reference/fdata.cen.md)
  : Center functional data
- [`deriv()`](https://sipemu.github.io/fdars-r/reference/deriv.md) :
  Compute functional derivative
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
- [`magnitudeshape()`](https://sipemu.github.io/fdars-r/reference/magnitudeshape.md)
  : Magnitude-Shape Outlier Detection for Functional Data

## Regression

- [`fregre.basis()`](https://sipemu.github.io/fdars-r/reference/fregre.basis.md)
  : Functional Basis Regression
- [`fregre.basis.cv()`](https://sipemu.github.io/fdars-r/reference/fregre.basis.cv.md)
  : Cross-Validation for Functional Basis Regression
- [`fregre.np()`](https://sipemu.github.io/fdars-r/reference/fregre.np.md)
  : Nonparametric Functional Regression
- [`fregre.np.cv()`](https://sipemu.github.io/fdars-r/reference/fregre.np.cv.md)
  : Cross-Validation for Nonparametric Functional Regression
- [`fregre.np.multi()`](https://sipemu.github.io/fdars-r/reference/fregre.np.multi.md)
  : Nonparametric Regression with Multiple Functional Predictors
- [`fregre.pc()`](https://sipemu.github.io/fdars-r/reference/fregre.pc.md)
  : Functional Regression
- [`fregre.pc.cv()`](https://sipemu.github.io/fdars-r/reference/fregre.pc.cv.md)
  : Cross-Validation for Functional PC Regression
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
- [`plot(`*`<cluster.kmeans>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.cluster.kmeans.md)
  : Plot Method for cluster.kmeans Objects
- [`plot(`*`<cluster.optim>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.cluster.optim.md)
  : Plot Method for cluster.optim Objects
- [`plot(`*`<fequiv.test>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.fequiv.test.md)
  : Plot method for fequiv.test
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
- [`plot(`*`<ssa_result>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.ssa_result.md)
  : Plot method for ssa_result objects
- [`plot(`*`<stl_result>`*`)`](https://sipemu.github.io/fdars-r/reference/plot.stl_result.md)
  : Plot method for stl_result objects

## Prediction

- [`predict(`*`<fregre.fd>`*`)`](https://sipemu.github.io/fdars-r/reference/predict.fregre.fd.md)
  : Predict Method for Functional Regression (fregre.fd)
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
- [`print(`*`<cluster.kmeans>`*`)`](https://sipemu.github.io/fdars-r/reference/print.cluster.kmeans.md)
  : Print Method for cluster.kmeans Objects
- [`print(`*`<cluster.optim>`*`)`](https://sipemu.github.io/fdars-r/reference/print.cluster.optim.md)
  : Print Method for cluster.optim Objects
- [`print(`*`<fbplot>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fbplot.md)
  : Print Method for fbplot Objects
- [`print(`*`<fregre.fd>`*`)`](https://sipemu.github.io/fdars-r/reference/print.fregre.fd.md)
  : Print method for fregre objects
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

## Internal (Rust Bindings)

Low-level Rust function wrappers. Use the R-level functions instead.

- [`alignment_align_to_target()`](https://sipemu.github.io/fdars-r/reference/alignment_align_to_target.md)
  : Align all curves to a target curve
- [`alignment_amplitude_dist()`](https://sipemu.github.io/fdars-r/reference/alignment_amplitude_dist.md)
  : Amplitude self-distance matrix
- [`alignment_compose_warps()`](https://sipemu.github.io/fdars-r/reference/alignment_compose_warps.md)
  : Compose two warping functions
- [`alignment_constrained()`](https://sipemu.github.io/fdars-r/reference/alignment_constrained.md)
  : Elastic alignment with landmark constraints
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
  : Full TSRVF transform
- [`alignment_warp_complexity()`](https://sipemu.github.io/fdars-r/reference/alignment_warp_complexity.md)
  : Compute warp complexity
- [`alignment_warp_smoothness()`](https://sipemu.github.io/fdars-r/reference/alignment_warp_smoothness.md)
  : Compute warp smoothness
- [`alignment_with_landmarks()`](https://sipemu.github.io/fdars-r/reference/alignment_with_landmarks.md)
  : Elastic alignment with automatic landmark detection
- [`landmark_detect()`](https://sipemu.github.io/fdars-r/reference/landmark_detect.md)
  : Detect landmarks in a single curve
- [`landmark_register_curves()`](https://sipemu.github.io/fdars-r/reference/landmark_register_curves.md)
  : Detect landmarks and register curves
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
- [`streaming_depth_batch()`](https://sipemu.github.io/fdars-r/reference/streaming_depth_batch.md)
  : Streaming depth: batch self-depth computation
- [`streaming_depth_one()`](https://sipemu.github.io/fdars-r/reference/streaming_depth_one.md)
  : Streaming depth: single curve against reference
- [`streaming_depth_vs_ref()`](https://sipemu.github.io/fdars-r/reference/streaming_depth_vs_ref.md)
  : Streaming depth: new data against reference
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
- [`depth_rp_2d()`](https://sipemu.github.io/fdars-r/reference/depth_rp_2d.md)
  : Random projection depth for 2D functional data (surfaces) Projects
  surfaces to scalars using random projections
- [`depth_rt_1d()`](https://sipemu.github.io/fdars-r/reference/depth_rt_1d.md)
  : Compute random Tukey depth
- [`depth_rt_2d()`](https://sipemu.github.io/fdars-r/reference/depth_rt_2d.md)
  : Random Tukey depth for 2D functional data (surfaces)
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
- [`fdata2basis()`](https://sipemu.github.io/fdars-r/reference/fdata2basis.md)
  : Convert Functional Data to Basis Coefficients
- [`fdata2basis_1d()`](https://sipemu.github.io/fdars-r/reference/fdata2basis_1d.md)
  : Convert functional data to basis coefficients type: 0 = bspline, 1 =
  fourier
- [`fdata2basis_2d()`](https://sipemu.github.io/fdars-r/reference/fdata2basis_2d.md)
  : Convert 2D Functional Data to Tensor Product Basis Coefficients
- [`fdata2basis_2d_raw()`](https://sipemu.github.io/fdars-r/reference/fdata2basis_2d_raw.md)
  : Project 2D functional data to tensor product basis coefficients (raw
  binding)
- [`fdata2basis_cv()`](https://sipemu.github.io/fdars-r/reference/fdata2basis_cv.md)
  : Cross-Validation for Basis Function Number Selection
- [`fdata2fd()`](https://sipemu.github.io/fdars-r/reference/fdata2fd.md)
  : Convert Functional Data to fd class
- [`fdata2pc()`](https://sipemu.github.io/fdars-r/reference/fdata2pc.md)
  : Convert Functional Data to Principal Component Scores
- [`fdata2pc_1d()`](https://sipemu.github.io/fdars-r/reference/fdata2pc_1d.md)
  : Perform functional PCA via SVD on centered data Returns: singular
  values, rotation matrix (loadings), scores, mean
- [`fdata2pls()`](https://sipemu.github.io/fdars-r/reference/fdata2pls.md)
  : Convert Functional Data to PLS Scores
- [`fdata2pls_1d()`](https://sipemu.github.io/fdars-r/reference/fdata2pls_1d.md)
  : Perform PLS via NIPALS algorithm Returns: weights, scores, loadings
- [`basis.aic()`](https://sipemu.github.io/fdars-r/reference/basis.aic.md)
  : AIC for Basis Representation
- [`basis.bic()`](https://sipemu.github.io/fdars-r/reference/basis.bic.md)
  : BIC for Basis Representation
- [`basis.gcv()`](https://sipemu.github.io/fdars-r/reference/basis.gcv.md)
  : GCV Score for Basis Representation
- [`basis2fdata()`](https://sipemu.github.io/fdars-r/reference/basis2fdata.md)
  : Basis Representation Functions for Functional Data
- [`basis2fdata_1d()`](https://sipemu.github.io/fdars-r/reference/basis2fdata_1d.md)
  : Reconstruct functional data from basis coefficients Returns data
  matrix (n x m)
- [`basis2fdata_2d()`](https://sipemu.github.io/fdars-r/reference/basis2fdata_2d.md)
  : Reconstruct 2D Functional Data from Tensor Product Basis
  Coefficients
- [`basis2fdata_2d_raw()`](https://sipemu.github.io/fdars-r/reference/basis2fdata_2d_raw.md)
  : Reconstruct 2D functional data from tensor product basis
  coefficients (raw binding)
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
- [`compute_adot()`](https://sipemu.github.io/fdars-r/reference/compute_adot.md)
  : Compute the Adot matrix (parallelized)
- [`pcvm_statistic()`](https://sipemu.github.io/fdars-r/reference/pcvm_statistic.md)
  : Compute the PCvM statistic
- [`rp_stat()`](https://sipemu.github.io/fdars-r/reference/rp_stat.md) :
  Compute random projection statistics (parallelized over projections)
- [`outliers_lrt()`](https://sipemu.github.io/fdars-r/reference/outliers_lrt.md)
  : LRT-based outlier detection Returns indices of detected outliers
- [`outliers_thres_lrt()`](https://sipemu.github.io/fdars-r/reference/outliers_thres_lrt.md)
  : Compute bootstrap threshold for LRT outlier detection Highly
  parallelized across bootstrap iterations
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
- [`kmeans_fd()`](https://sipemu.github.io/fdars-r/reference/kmeans_fd.md)
  : Functional k-means clustering
- [`fuzzycmeans_fd()`](https://sipemu.github.io/fdars-r/reference/fuzzycmeans_fd.md)
  : Fuzzy C-Means clustering for functional data m_fuzz is the fuzziness
  parameter (typically 2)
- [`register_shift_1d()`](https://sipemu.github.io/fdars-r/reference/register_shift_1d.md)
  : Shift registration: find optimal horizontal shift for each curve to
  align with a target (usually the mean)
- [`int_simpson()`](https://sipemu.github.io/fdars-r/reference/int_simpson.md)
  : Simpson's rule integration for functional data Integrates each curve
  over the domain
- [`inprod_fdata()`](https://sipemu.github.io/fdars-r/reference/inprod_fdata.md)
  : Inner product of two functional data objects \<f, g\> =
  integral(f(t) \* g(t) dt)
- [`knn_gcv()`](https://sipemu.github.io/fdars-r/reference/knn_gcv.md) :
  k-NN with Global Cross-Validation Finds a single optimal k for all
  observations
- [`knn_lcv()`](https://sipemu.github.io/fdars-r/reference/knn_lcv.md) :
  k-NN with Local Cross-Validation Finds an optimal k for each
  observation
- [`knn_predict()`](https://sipemu.github.io/fdars-r/reference/knn_predict.md)
  : Kernel prediction with fixed bandwidth for prediction on new data
- [`silhouette_score()`](https://sipemu.github.io/fdars-r/reference/silhouette_score.md)
  : Compute silhouette score for clustering Returns the mean silhouette
  coefficient across all samples
- [`calinski_harabasz()`](https://sipemu.github.io/fdars-r/reference/calinski_harabasz.md)
  : Compute Calinski-Harabasz index (variance ratio criterion) Higher
  values indicate better defined clusters
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
- [`eigenfunctions_1d()`](https://sipemu.github.io/fdars-r/reference/eigenfunctions_1d.md)
  : Compute eigenfunction basis values efun_type: 0 = Fourier, 1 = Poly,
  2 = PolyHigh, 3 = Wiener
- [`eigenvalues_1d()`](https://sipemu.github.io/fdars-r/reference/eigenvalues_1d.md)
  : Generate eigenvalue sequence eval_type: 0 = linear, 1 = exponential,
  2 = wiener
- [`sim_kl_1d()`](https://sipemu.github.io/fdars-r/reference/sim_kl_1d.md)
  : Simulate functional data via Karhunen-Loève expansion
- [`add_error_curve_1d()`](https://sipemu.github.io/fdars-r/reference/add_error_curve_1d.md)
  : Add curve-level Gaussian noise to functional data
- [`add_error_pointwise_1d()`](https://sipemu.github.io/fdars-r/reference/add_error_pointwise_1d.md)
  : Add pointwise Gaussian noise to functional data
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
- [`select_basis_auto()`](https://sipemu.github.io/fdars-r/reference/select_basis_auto.md)
  : Automatic basis selection for each curve individually.
- [`pspline_fit_1d()`](https://sipemu.github.io/fdars-r/reference/pspline_fit_1d.md)
  : P-spline fitting: returns coefficients, fitted values, and
  diagnostics
- [`pspline_fit_2d()`](https://sipemu.github.io/fdars-r/reference/pspline_fit_2d.md)
  : 2D P-spline fitting with anisotropic penalties
- [`geometric_median_1d()`](https://sipemu.github.io/fdars-r/reference/geometric_median_1d.md)
  : Compute the geometric median (L1 median) of functional data using
  Weiszfeld's algorithm The geometric median minimizes sum of L2
  distances to all curves
- [`geometric_median_2d()`](https://sipemu.github.io/fdars-r/reference/geometric_median_2d.md)
  : Compute the geometric median (L1 median) of 2D functional data using
  Weiszfeld's algorithm Data is stored as n x (m1\*m2) matrix where each
  row is a flattened surface
