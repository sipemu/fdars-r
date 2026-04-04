# CUSUM Monitoring for Functional Data

Applies Cumulative Sum (CUSUM) monitoring to sequential functional data.
CUSUM accumulates evidence of a mean shift over time, providing faster
detection of small persistent shifts than Shewhart-type charts.

## Usage

``` r
spm.cusum(
  chart,
  newdata,
  k = 0.5,
  h = 5,
  ncomp = NULL,
  alpha = 0.05,
  restart = FALSE
)
```

## Arguments

- chart:

  An object of class `spm.chart` from
  [`spm.phase1`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md).

- newdata:

  An object of class `fdata` with sequential observations.

- k:

  CUSUM allowance (slack) parameter (default 0.5). Controls the size of
  shift the CUSUM is tuned to detect.

- h:

  CUSUM decision interval (default 5.0). Larger values give fewer false
  alarms but slower detection.

- ncomp:

  Number of components for CUSUM (default: same as chart).

- alpha:

  Significance level for SPE limit (default 0.05).

- restart:

  Logical; if `TRUE`, the CUSUM resets to zero after each alarm (default
  `FALSE`).

## Value

An object of class `spm.cusum` with components:

- cusum.statistic:

  Numeric vector of CUSUM values

- ucl:

  Decision interval (h)

- alarm:

  Logical vector: TRUE where CUSUM exceeds h

- scores:

  FPC score matrix

- spe:

  SPE values

- spe.limit:

  SPE control limit

- spe.alarm:

  Logical: TRUE where SPE exceeds limit

## See also

[`spm.monitor`](https://sipemu.github.io/fdars-r/reference/spm.monitor.md)
for Shewhart-type monitoring,
[`spm.mewma`](https://sipemu.github.io/fdars-r/reference/spm.mewma.md)
for MEWMA monitoring

## Examples

``` r
# \donttest{
set.seed(1)
n <- 50; m <- 30
argvals <- seq(0, 1, length.out = m)
X <- matrix(rnorm(n * m), n, m)
fd <- fdata(X, argvals = argvals)
chart <- spm.phase1(fd, ncomp = 3)

# Monitor with a small persistent shift
X_new <- matrix(rnorm(30 * m) + 0.5, 30, m)
fd_new <- fdata(X_new, argvals = argvals)
cusum <- spm.cusum(chart, fd_new, k = 0.5, h = 5.0)
cusum
#> SPM CUSUM Monitoring
#>   Observations: 30 
#>   k: 0.5 
#>   h (UCL): 5 
#>   Restart: FALSE 
#>   CUSUM alarms: 22 of 30 (73.3%) 
#>   SPE alarms: 2 of 30 (6.67%) 
# }
```
