# Iterative Phase I Chart Construction

Builds a Phase I control chart iteratively by removing observations that
exceed control limits, then re-fitting. This produces cleaner reference
limits when the training data may contain outliers.

## Usage

``` r
spm.phase1.iterative(
  fdataobj,
  ncomp = 5,
  alpha = 0.05,
  tuning.fraction = 0.5,
  seed = 42,
  max.iterations = 10,
  max.removal.fraction = 0.3
)
```

## Arguments

- fdataobj:

  An object of class `fdata` with in-control data.

- ncomp:

  Number of principal components (default 5).

- alpha:

  Significance level (default 0.05).

- tuning.fraction:

  Fraction for FPCA tuning (default 0.5).

- seed:

  Random seed (default 42).

- max.iterations:

  Maximum number of iterative removal rounds (default 10).

- max.removal.fraction:

  Maximum fraction of data that may be removed (default 0.3).

## Value

An object of class `spm.chart` (compatible with `spm.monitor`,
`spm.cusum`, etc.) with additional fields:

- n.iterations:

  Number of iterations performed

- removed.indices:

  Integer vector of removed observation indices (1-indexed)

- n.remaining:

  Number of observations remaining

- removal.history:

  List of integer vectors; removed indices per iteration

- removal.rates:

  Numeric vector of removal rates per iteration

## See also

[`spm.phase1`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md)
for non-iterative Phase I

## Examples

``` r
# \donttest{
set.seed(1)
n <- 60; m <- 30
argvals <- seq(0, 1, length.out = m)
# Add some outliers
X <- matrix(rnorm(n * m), n, m)
X[1:3, ] <- X[1:3, ] + 5
fd <- fdata(X, argvals = argvals)

chart <- spm.phase1.iterative(fd, ncomp = 3, max.iterations = 5)
chart
#> SPM Control Chart (Phase I)
#>   Components: 3 
#>   Alpha: 0.05 
#>   T2 UCL: 7.815 
#>   SPE UCL: 1.366 
#>   Observations: 60 
#>   Grid points: 30 
#>   Eigenvalues: 3.96, 3.21, 2.99 
chart$removed.indices
#> [1]  1  2  3 46 34 15 22 48 35
# }
```
