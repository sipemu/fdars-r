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
# \donttest{
set.seed(1)
n <- 10
n_proj <- 3
proj <- unlist(replicate(n_proj, sample(n), simplify = FALSE))
rp.stat(proj, rnorm(n), n.proj = n_proj)
#> $cvm
#> [1] 0.05562997 0.16858091 0.19217356
#> 
#> $ks
#> [1] 0.3426421 0.6733942 0.7741449
#> 
# }
```
