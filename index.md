[![R-CMD-check](https://github.com/sipemu/fdars-r/actions/workflows/r-cmd-check.yml/badge.svg)](https://github.com/sipemu/fdars-r/actions/workflows/r-cmd-check.yml)
[![CRAN
status](https://www.r-pkg.org/badges/version/fdars)](https://CRAN.R-project.org/package=fdars)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/fdars)](https://cran.r-project.org/package=fdars)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![codecov](https://codecov.io/gh/sipemu/fdars-r/graph/badge.svg)](https://app.codecov.io/gh/sipemu/fdars-r)

## Functional Data Analysis in Rust

**fdars** is a comprehensive R toolkit for functional data analysis
powered by a high-performance Rust backend. Treat entire curves,
spectra, and trajectories as single observations — then smooth, align,
decompose, and analyze them with a consistent, pipe-friendly API.

[Get started](https://sipemu.github.io/fdars-r/articles/introduction.md)
[Reference](https://sipemu.github.io/fdars-r/reference/index.md)

------------------------------------------------------------------------

## Learn

Get started, preprocess, and simulate functional data.

[](https://sipemu.github.io/fdars-r/articles/introduction.md)

![Introduction](reference/figures/card-introduction.svg)

### Introduction to fdars

Get started with functional data objects and basic operations

[](https://sipemu.github.io/fdars-r/articles/custom-plotting.md)

![Custom Plotting](reference/figures/card-custom-plotting.svg)

### Custom Plotting

Customize curve visualizations with ggplot2

[](https://sipemu.github.io/fdars-r/articles/intro-to-smoothing.md)

![Smoothing](reference/figures/card-intro-to-smoothing.svg)

### Introduction to Smoothing

Smooth noisy curves using kernel and spline methods

[](https://sipemu.github.io/fdars-r/articles/working-with-derivatives.md)

![Derivatives](reference/figures/card-working-with-derivatives.svg)

### Working with Derivatives

Compute and analyze functional derivatives

[](https://sipemu.github.io/fdars-r/articles/simulation-toolbox.md)

![Simulation](reference/figures/card-simulation-toolbox.svg)

### Simulation Toolbox

Generate synthetic functional data

[](https://sipemu.github.io/fdars-r/articles/irregular-sampling.md)

![Irregular Sampling](reference/figures/card-irregular-sampling.svg)

### Irregular Sampling

Work with irregularly sampled and sparse data

## Represent

Decompose, transform, rank, and measure functional data.

[](https://sipemu.github.io/fdars-r/articles/fpca.md)

![FPCA](reference/figures/card-fpca.svg)

### Functional PCA

Extract dominant modes of variation

[](https://sipemu.github.io/fdars-r/articles/basis-representation.md)

![Basis Representation](reference/figures/card-basis-representation.svg)

### Basis Representation

Expand curves in B-spline and Fourier bases

[](https://sipemu.github.io/fdars-r/articles/andrews-transformation.md)

![Andrews Curves](reference/figures/card-andrews-transformation.svg)

### Andrews Curves

Transform multivariate data into functional curves

[](https://sipemu.github.io/fdars-r/articles/depth-functions.md)

![Depth Functions](reference/figures/card-depth-functions.svg)

### Depth Functions

Rank curves from center outward using statistical depth

[](https://sipemu.github.io/fdars-r/articles/streaming-depth.md)

![Streaming Depth](reference/figures/card-streaming-depth.svg)

### Streaming Depth

Monitor depth in real-time as new curves arrive

[](https://sipemu.github.io/fdars-r/articles/distance-metrics.md)

![Distance Metrics](reference/figures/card-distance-metrics.svg)

### Distance Metrics

Measure similarity with L^(p), DTW, and elastic distances

## Align

Register and align curves to remove phase variability.

[](https://sipemu.github.io/fdars-r/articles/elastic-alignment.md)

![Elastic Alignment](reference/figures/card-elastic-alignment.svg)

### Elastic Alignment

Remove phase variability via SRSF

[](https://sipemu.github.io/fdars-r/articles/landmark-registration.md)

![Landmark
Registration](reference/figures/card-landmark-registration.svg)

### Landmark Registration

Align curves by matching peaks and valleys

[](https://sipemu.github.io/fdars-r/articles/tsrvf.md)

![TSRVF](reference/figures/card-tsrvf.svg)

### TSRVF (Tangent Space)

Project aligned curves into a linear tangent space

[](https://sipemu.github.io/fdars-r/articles/alignment-comparison.md)

![Comparing Methods](reference/figures/card-alignment-comparison.svg)

### Comparing Methods

Compare alignment methods side-by-side

## Analyze

Infer, classify, predict, and model functional data.

[](https://sipemu.github.io/fdars-r/articles/tolerance-bands.md)

![Tolerance Bands](reference/figures/card-tolerance-bands.svg)

### Tolerance Bands

Construct tolerance and confidence bands

[](https://sipemu.github.io/fdars-r/articles/equivalence-testing.md)

![Equivalence Testing](reference/figures/card-equivalence-testing.svg)

### Equivalence Testing

Test functional equivalence between groups

[](https://sipemu.github.io/fdars-r/articles/clustering.md)

![Clustering](reference/figures/card-clustering.svg)

### Clustering

Group similar curves with k-means and fuzzy c-means

[](https://sipemu.github.io/fdars-r/articles/outlier-detection.md)

![Outlier Detection](reference/figures/card-outlier-detection.svg)

### Outlier Detection

Identify anomalous curves using depth

[](https://sipemu.github.io/fdars-r/articles/regression.md)

![Regression](reference/figures/card-regression.svg)

### Regression

Predict scalar outcomes from functional predictors

[](https://sipemu.github.io/fdars-r/articles/seasonal-analysis.md)

![Seasonal Analysis](reference/figures/card-seasonal-analysis.svg)

### Seasonal Analysis

Detect and decompose seasonal patterns

[](https://sipemu.github.io/fdars-r/articles/covariance-functions.md)

![Covariance Functions](reference/figures/card-covariance-functions.svg)

### Covariance Functions

Build GP models with composable kernels
