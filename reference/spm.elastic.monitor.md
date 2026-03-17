# Elastic SPM Phase II Monitoring

Monitors new functional data against an elastic SPM chart. New curves
are aligned to the Karcher mean, then amplitude (and optionally phase)
are monitored separately.

## Usage

``` r
spm.elastic.monitor(chart, newdata)
```

## Arguments

- chart:

  An object of class `spm.elastic.chart` from
  [`spm.elastic.phase1`](https://sipemu.github.io/fdars-r/reference/spm.elastic.phase1.md).

- newdata:

  An object of class `fdata` with new observations.

## Value

A list with components:

- amplitude:

  List with t2, spe, t2.alarm, spe.alarm, scores

- phase:

  List with t2, spe, t2.alarm, spe.alarm, scores (if
  `monitor.phase = TRUE`), or NULL

- aligned.data:

  Matrix of aligned curves

- warping.functions:

  Matrix of estimated warping functions

## See also

[`spm.elastic.phase1`](https://sipemu.github.io/fdars-r/reference/spm.elastic.phase1.md)
for building the chart

## Examples

``` r
# \donttest{
set.seed(1)
n <- 40; m <- 30
argvals <- seq(0, 1, length.out = m)
X <- matrix(rnorm(n * m), n, m)
fd <- fdata(X, argvals = argvals)
chart <- spm.elastic.phase1(fd, ncomp = 3, warp.ncomp = 2)

X_new <- matrix(rnorm(10 * m) + 1, 10, m)
fd_new <- fdata(X_new, argvals = argvals)
mon <- spm.elastic.monitor(chart, fd_new)
mon$amplitude$t2.alarm
#>  [1] FALSE  TRUE FALSE FALSE  TRUE  TRUE FALSE FALSE  TRUE  TRUE
# }
```
