# Elastic Changepoint Detection

Detect structural breaks in functional time series using elastic
methods. Elastic Changepoint Detection

## Usage

``` r
elastic.changepoint(
  fdataobj,
  type = c("amplitude", "phase", "fpca"),
  pca.method = c("vertical", "horizontal", "joint"),
  ncomp = 3,
  lambda = 0,
  max.iter = 20,
  n.mc = 1000,
  cov.kernel = c("bartlett", "parzen", "flattop", "simple"),
  cov.bandwidth = NULL,
  seed = NULL
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'. Curves should be in temporal order (rows =
  time points in the series).

- type:

  Type of changepoint test: "amplitude" (shape changes), "phase" (timing
  changes), or "fpca" (combined via elastic FPCA).

- pca.method:

  PCA method for type = "fpca": "vertical", "horizontal", or "joint".

- ncomp:

  Number of FPC components for type = "fpca" (default 3).

- lambda:

  Regularization for alignment (default 0).

- max.iter:

  Maximum alignment iterations (default 20).

- n.mc:

  Number of Monte Carlo permutations for p-value (default 1000).

- cov.kernel:

  Long-run covariance kernel: "bartlett", "parzen", "flattop", or
  "simple".

- cov.bandwidth:

  Bandwidth for the covariance kernel. If NULL, automatically selected.

- seed:

  Random seed for reproducibility.

## Value

An object of class 'elastic.changepoint' with components:

- changepoint:

  Estimated changepoint index (1-based)

- test.statistic:

  Value of the CUSUM test statistic

- p.value:

  Permutation-based p-value

- cusum.values:

  Full CUSUM process values

- type:

  Type of test performed

- n:

  Number of curves in the series

## Details

Tests for a changepoint in a functional time series using amplitude,
phase, or FPCA-based test statistics with permutation-based p-values.
