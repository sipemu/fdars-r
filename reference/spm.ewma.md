# EWMA-Based SPM Monitoring

Applies Exponentially Weighted Moving Average (EWMA) smoothing to FPC
scores before computing monitoring statistics. This increases
sensitivity to small persistent shifts in the process.

## Usage

``` r
spm.ewma(chart, newdata, lambda = 0.2, ncomp = NULL, alpha = 0.05)
```

## Arguments

- chart:

  An object of class `spm.chart` from
  [`spm.phase1`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md).

- newdata:

  An object of class `fdata` with sequential observations (rows in time
  order).

- lambda:

  EWMA smoothing parameter in (0, 1\] (default 0.2). Smaller values give
  more smoothing; lambda = 1 gives raw scores.

- ncomp:

  Number of components for EWMA (default: same as chart).

- alpha:

  Significance level for EWMA control limit (default 0.05).

## Value

A list with components:

- smoothed.scores:

  EWMA-smoothed score matrix

- t2:

  T-squared values on smoothed scores

- spe:

  SPE values (unsmoothed)

- t2.limit:

  T-squared control limit for EWMA

- spe.limit:

  SPE control limit

- t2.alarm:

  Logical: TRUE where T-squared exceeds limit

- spe.alarm:

  Logical: TRUE where SPE exceeds limit

## See also

[`spm.phase1`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md)
for building the chart,
[`spm.monitor`](https://sipemu.github.io/fdars-r/reference/spm.monitor.md)
for standard (non-EWMA) monitoring

## Examples

``` r
# \donttest{
set.seed(1)
n <- 50; m <- 30
argvals <- seq(0, 1, length.out = m)
X <- matrix(rnorm(n * m), n, m)
fd <- fdata(X, argvals = argvals)
chart <- spm.phase1(fd, ncomp = 3)

# Sequential monitoring with small shift
X_seq <- matrix(rnorm(20 * m) + 0.5, 20, m)
fd_seq <- fdata(X_seq, argvals = argvals)
ewma_result <- spm.ewma(chart, fd_seq, lambda = 0.2)
ewma_result$t2.alarm
#>  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#> [13] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# }
```
