# Gaussian Mixture Model Clustering for Functional Data

Performs model-based clustering using Gaussian mixture models on
functional data. Curves are projected onto a B-spline or Fourier basis
before fitting. Automatic K selection via BIC or ICL.

## Usage

``` r
cluster.gmm(
  fdataobj,
  k.range = 2:6,
  covariates = NULL,
  nbasis = 10,
  basis.type = c("bspline", "fourier"),
  cov.type = c("full", "diagonal"),
  cov.weight = 1,
  max.iter = 100,
  tol = 1e-04,
  n.init = 5,
  seed = NULL,
  criterion = c("bic", "icl")
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- k.range:

  Range of number of clusters to try (default 2:6).

- covariates:

  Optional matrix of scalar covariates to include.

- nbasis:

  Number of basis functions for projection (default 10).

- basis.type:

  Basis type: "bspline" (default) or "fourier".

- cov.type:

  Covariance structure: "full" (default) or "diagonal".

- cov.weight:

  Weight for covariates vs basis coefficients (default 1).

- max.iter:

  Maximum EM iterations (default 100).

- tol:

  Convergence tolerance (default 1e-4).

- n.init:

  Number of random initializations (default 5).

- seed:

  Optional random seed.

- criterion:

  Model selection criterion: "bic" (default) or "icl".

## Value

An object of class 'cluster.gmm' with components:

- cluster:

  Integer vector of cluster assignments (1-indexed)

- membership:

  Matrix of posterior membership probabilities

- means:

  Component means matrix

- weights:

  Mixing proportions

- bic:

  BIC of selected model

- icl:

  ICL of selected model

- k:

  Number of clusters selected

- converged:

  Whether EM converged

- bic.values:

  BIC values for each K tried

## See also

[`cluster.kmeans`](https://sipemu.github.io/fdars-r/reference/cluster.kmeans.md)
for hard clustering,
[`cluster.fcm`](https://sipemu.github.io/fdars-r/reference/cluster.fcm.md)
for fuzzy clustering
