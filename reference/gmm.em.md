# Gaussian Mixture Model EM on Feature Matrix

Fits a GMM directly to a feature matrix using the EM algorithm. Unlike
[`cluster.gmm`](https://sipemu.github.io/fdars-r/reference/cluster.gmm.md)
which operates on functional data, this function works on arbitrary
numeric feature matrices.

## Usage

``` r
gmm.em(
  features,
  k,
  cov.type = c("full", "diagonal", "spherical"),
  max.iter = 100,
  tol = 1e-06,
  seed = 42
)
```

## Arguments

- features:

  Numeric matrix of features (n x d).

- k:

  Number of mixture components.

- cov.type:

  Covariance type: "full", "diagonal", or "spherical".

- max.iter:

  Maximum EM iterations (default 100).

- tol:

  Convergence tolerance (default 1e-6).

- seed:

  Random seed.

## Value

A list with components:

- `cluster` — Integer vector of hard cluster assignments.

- `membership` — Posterior membership probability matrix (n x k).

- `means` — List of component mean vectors.

- `weights` — Mixing proportions.

- `log.likelihood` — Final log-likelihood.

- `bic`, `icl` — Model selection criteria.

- `iterations` — Number of EM iterations.

## Examples

``` r
# \donttest{
X <- rbind(matrix(rnorm(100), 50, 2), matrix(rnorm(100, 3), 50, 2))
fit <- gmm.em(X, k = 2)
table(fit$cluster)
#> 
#>  2  3 
#> 50 50 
# }
```
