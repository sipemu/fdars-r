# Monitor a Batch of Partially-Observed Curves

Monitors multiple partially-observed curves in a single call. Each curve
may have a different number of observed grid points.

## Usage

``` r
spm.monitor.partial.batch(
  chart,
  partial.data.list,
  n.observed.list,
  completion = c("conditional", "projection", "zero_pad"),
  ncomp = NULL,
  alpha = 0.05
)
```

## Arguments

- chart:

  An object of class `spm.chart` from
  [`spm.phase1`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md).

- partial.data.list:

  A list of numeric vectors, each containing the observed values for one
  curve.

- n.observed.list:

  An integer vector giving the number of observed grid points for each
  curve.

- completion:

  Completion method: `"conditional"`, `"projection"`, or `"zero_pad"`.

- ncomp:

  Number of components (default: same as chart).

- alpha:

  Significance level (default 0.05).

## Value

A list of results, one per curve. Each element is a list with:

- scores:

  Estimated FPC scores

- t2:

  T-squared statistic

- t2.alarm:

  Logical

- domain.fraction:

  Fraction of domain observed

- completed.curve:

  The completed curve, or NULL

## See also

[`spm.monitor.partial`](https://sipemu.github.io/fdars-r/reference/spm.monitor.partial.md)
for a single curve

## Examples

``` r
# \donttest{
set.seed(1)
n <- 50; m <- 30
argvals <- seq(0, 1, length.out = m)
X <- matrix(rnorm(n * m), n, m)
fd <- fdata(X, argvals = argvals)
chart <- spm.phase1(fd, ncomp = 3)

# Three curves with different observed lengths
curves <- list(rnorm(30), rnorm(30), rnorm(30))
n_obs <- c(10L, 20L, 25L)
results <- spm.monitor.partial.batch(chart, curves, n_obs)
sapply(results, function(r) r$t2)
#> [1] 3.300021 1.962663 1.282167
# }
```
