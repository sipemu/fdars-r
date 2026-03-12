# Random Projection Statistic

Computes Cramer-von Mises and Kolmogorov-Smirnov statistics for random
projection goodness-of-fit testing.

## Usage

``` r
rp.stat(proj.x.ord, residuals, n.proj)
```

## Arguments

- proj.x.ord:

  Integer vector of projection orderings.

- residuals:

  Numeric vector of residuals.

- n.proj:

  Number of random projections.

## Value

A list with `cvm` and `ks` vectors of test statistics.

## Examples

``` r
rp.stat(sample(10), rnorm(10), n.proj = 3)
#> Error in rp_stat(as.integer(proj.x.ord), as.double(residuals), as.integer(n.proj)): User function panicked: rp_stat
```
