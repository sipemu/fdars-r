---
output: github_document
---

<!-- index.md is the pkgdown home page; README.md remains for GitHub -->

[![R-CMD-check](https://github.com/sipemu/fdars-r/actions/workflows/r-cmd-check.yml/badge.svg)](https://github.com/sipemu/fdars-r/actions/workflows/r-cmd-check.yml)
[![CRAN status](https://www.r-pkg.org/badges/version/fdars)](https://CRAN.R-project.org/package=fdars)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/fdars)](https://cran.r-project.org/package=fdars)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![codecov](https://codecov.io/gh/sipemu/fdars-r/graph/badge.svg)](https://app.codecov.io/gh/sipemu/fdars-r)

## Functional Data Analysis in Rust

**fdars** is a comprehensive R toolkit for functional data analysis powered by a high-performance Rust backend. Treat entire curves, spectra, and trajectories as single observations — then smooth, align, decompose, and analyze them with a consistent, pipe-friendly API.

---

```{=html}
<!-- ===== LEARN ===== -->
<h2 class="fdars-section-heading fdars-learn">Learn</h2>
<p class="fdars-section-desc">Get started, preprocess, and simulate functional data.</p>
<div class="fdars-gallery">

  <a class="fdars-gallery-item" href="articles/introduction.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-introduction.svg" alt="Introduction"></div>
    <div class="fdars-gallery-title">Introduction to fdars</div>
  </a>

  <a class="fdars-gallery-item" href="articles/custom-plotting.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-custom-plotting.svg" alt="Custom Plotting"></div>
    <div class="fdars-gallery-title">Custom Plotting</div>
  </a>

  <a class="fdars-gallery-item" href="articles/intro-to-smoothing.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-intro-to-smoothing.svg" alt="Smoothing"></div>
    <div class="fdars-gallery-title">Introduction to Smoothing</div>
  </a>

  <a class="fdars-gallery-item" href="articles/working-with-derivatives.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-working-with-derivatives.svg" alt="Derivatives"></div>
    <div class="fdars-gallery-title">Working with Derivatives</div>
  </a>

  <a class="fdars-gallery-item" href="articles/simulation-toolbox.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-simulation-toolbox.svg" alt="Simulation"></div>
    <div class="fdars-gallery-title">Simulation Toolbox</div>
  </a>

  <a class="fdars-gallery-item" href="articles/irregular-sampling.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-irregular-sampling.svg" alt="Irregular Sampling"></div>
    <div class="fdars-gallery-title">Irregular Sampling</div>
  </a>

</div>

<!-- ===== REPRESENT ===== -->
<h2 class="fdars-section-heading fdars-represent">Represent</h2>
<p class="fdars-section-desc">Decompose, transform, rank, and measure functional data.</p>
<div class="fdars-gallery">

  <a class="fdars-gallery-item" href="articles/fpca.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-fpca.svg" alt="FPCA"></div>
    <div class="fdars-gallery-title">Functional PCA</div>
  </a>

  <a class="fdars-gallery-item" href="articles/basis-representation.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-basis-representation.svg" alt="Basis Representation"></div>
    <div class="fdars-gallery-title">Basis Representation</div>
  </a>

  <a class="fdars-gallery-item" href="articles/andrews-transformation.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-andrews-transformation.svg" alt="Andrews Curves"></div>
    <div class="fdars-gallery-title">Andrews Curves</div>
  </a>

  <a class="fdars-gallery-item" href="articles/depth-functions.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-depth-functions.svg" alt="Depth Functions"></div>
    <div class="fdars-gallery-title">Depth Functions</div>
  </a>

  <a class="fdars-gallery-item" href="articles/streaming-depth.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-streaming-depth.svg" alt="Streaming Depth"></div>
    <div class="fdars-gallery-title">Streaming Depth</div>
  </a>

  <a class="fdars-gallery-item" href="articles/distance-metrics.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-distance-metrics.svg" alt="Distance Metrics"></div>
    <div class="fdars-gallery-title">Distance Metrics</div>
  </a>

  <a class="fdars-gallery-item" href="articles/elastic-fpca.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-elastic-fpca.svg" alt="Elastic FPCA"></div>
    <div class="fdars-gallery-title">Elastic FPCA</div>
  </a>

</div>

<!-- ===== ALIGN ===== -->
<h2 class="fdars-section-heading fdars-align">Align</h2>
<p class="fdars-section-desc">Register and align curves to remove phase variability.</p>
<div class="fdars-gallery">

  <a class="fdars-gallery-item" href="articles/elastic-alignment.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-elastic-alignment.svg" alt="Elastic Alignment"></div>
    <div class="fdars-gallery-title">Elastic Alignment</div>
  </a>

  <a class="fdars-gallery-item" href="articles/landmark-registration.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-landmark-registration.svg" alt="Landmark Registration"></div>
    <div class="fdars-gallery-title">Landmark Registration</div>
  </a>

  <a class="fdars-gallery-item" href="articles/tsrvf.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-tsrvf.svg" alt="TSRVF"></div>
    <div class="fdars-gallery-title">TSRVF (Tangent Space)</div>
  </a>

  <a class="fdars-gallery-item" href="articles/alignment-comparison.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-alignment-comparison.svg" alt="Comparing Methods"></div>
    <div class="fdars-gallery-title">Comparing Methods</div>
  </a>

</div>

<!-- ===== REGRESSION ===== -->
<h2 class="fdars-section-heading fdars-regression">Regression</h2>
<p class="fdars-section-desc">Regression, classification, and mixed models for functional data.</p>
<div class="fdars-gallery">

  <a class="fdars-gallery-item" href="articles/scalar-on-function.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-scalar-on-function.svg" alt="Scalar-on-Function Regression"></div>
    <div class="fdars-gallery-title">Scalar-on-Function Regression</div>
  </a>

  <a class="fdars-gallery-item" href="articles/function-on-scalar.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-function-on-scalar.svg" alt="Function-on-Scalar Regression"></div>
    <div class="fdars-gallery-title">Function-on-Scalar Regression</div>
  </a>

  <a class="fdars-gallery-item" href="articles/functional-classification.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-functional-classification.svg" alt="Classification"></div>
    <div class="fdars-gallery-title">Classification</div>
  </a>

  <a class="fdars-gallery-item" href="articles/functional-mixed-models.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-functional-mixed-models.svg" alt="Mixed Models"></div>
    <div class="fdars-gallery-title">Mixed Models</div>
  </a>

  <a class="fdars-gallery-item" href="articles/cross-validation.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-cross-validation.svg" alt="Cross-Validation"></div>
    <div class="fdars-gallery-title">Cross-Validation</div>
  </a>

  <a class="fdars-gallery-item" href="articles/elastic-regression.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-elastic-regression.svg" alt="Elastic Regression"></div>
    <div class="fdars-gallery-title">Elastic Regression</div>
  </a>

  <a class="fdars-gallery-item" href="articles/explainability.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-explainability.svg" alt="Model Explainability"></div>
    <div class="fdars-gallery-title">Model Explainability</div>
  </a>

  <a class="fdars-gallery-item" href="articles/regression-diagnostics.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-regression-diagnostics.svg" alt="Regression Diagnostics"></div>
    <div class="fdars-gallery-title">Regression Diagnostics</div>
  </a>

  <a class="fdars-gallery-item" href="articles/uncertainty-quantification.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-uncertainty-quantification.svg" alt="Uncertainty Quantification"></div>
    <div class="fdars-gallery-title">Uncertainty Quantification</div>
  </a>

  <a class="fdars-gallery-item" href="articles/smooth-basis.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-smooth-basis.svg" alt="Penalized Basis Smoothing"></div>
    <div class="fdars-gallery-title">Penalized Basis Smoothing</div>
  </a>

</div>

<!-- ===== ANALYZE ===== -->
<h2 class="fdars-section-heading fdars-analyze">Analyze</h2>
<p class="fdars-section-desc">Infer, cluster, detect outliers, and test functional data.</p>
<div class="fdars-gallery">

  <a class="fdars-gallery-item" href="articles/tolerance-bands.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-tolerance-bands.svg" alt="Tolerance Bands"></div>
    <div class="fdars-gallery-title">Tolerance Bands</div>
  </a>

  <a class="fdars-gallery-item" href="articles/equivalence-testing.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-equivalence-testing.svg" alt="Equivalence Testing"></div>
    <div class="fdars-gallery-title">Equivalence Testing</div>
  </a>

  <a class="fdars-gallery-item" href="articles/clustering.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-clustering.svg" alt="Clustering"></div>
    <div class="fdars-gallery-title">Clustering</div>
  </a>

  <a class="fdars-gallery-item" href="articles/gmm-clustering.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-gmm-clustering.svg" alt="GMM Clustering"></div>
    <div class="fdars-gallery-title">GMM Clustering</div>
  </a>

  <a class="fdars-gallery-item" href="articles/outlier-detection.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-outlier-detection.svg" alt="Outlier Detection"></div>
    <div class="fdars-gallery-title">Outlier Detection</div>
  </a>

  <a class="fdars-gallery-item" href="articles/seasonal-analysis.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-seasonal-analysis.svg" alt="Seasonal Analysis"></div>
    <div class="fdars-gallery-title">Seasonal Analysis</div>
  </a>

  <a class="fdars-gallery-item" href="articles/covariance-functions.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-covariance-functions.svg" alt="Covariance Functions"></div>
    <div class="fdars-gallery-title">Covariance Functions</div>
  </a>

  <a class="fdars-gallery-item" href="articles/elastic-changepoint.html">
    <div class="fdars-gallery-thumb"><img src="reference/figures/card-elastic-changepoint.svg" alt="Elastic Changepoint"></div>
    <div class="fdars-gallery-title">Elastic Changepoint</div>
  </a>

</div>
```
