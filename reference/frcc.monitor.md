# Monitor New Data Against FRCC Chart (Phase II)

Predicts the functional response from the FOSR model, computes
residuals, and monitors them against the established FRCC chart.

## Usage

``` r
frcc.monitor(chart, newdata, new.predictors)
```

## Arguments

- chart:

  An object of class `frcc.chart` from
  [`frcc.phase1`](https://sipemu.github.io/fdars-r/reference/frcc.phase1.md).

- newdata:

  An object of class `fdata` with new functional response.

- new.predictors:

  A matrix of new scalar predictors (n_new x p).

## Value

An object of class `spm.monitor` with components:

- t2:

  T-squared values

- spe:

  SPE values

- t2.alarm:

  Logical: TRUE where T-squared exceeds UCL

- spe.alarm:

  Logical: TRUE where SPE exceeds UCL

- residual.scores:

  Residual FPC scores

- t2.ucl:

  T-squared control limit

- spe.ucl:

  SPE control limit

## See also

[`frcc.phase1`](https://sipemu.github.io/fdars-r/reference/frcc.phase1.md)
for building the chart

## Examples

``` r
# \donttest{
set.seed(1)
n <- 60; m <- 20
argvals <- seq(0, 1, length.out = m)
X_pred <- cbind(rnorm(n), rnorm(n))
Y <- matrix(rnorm(n * m), n, m)
fd <- fdata(Y, argvals = argvals)
chart <- frcc.phase1(fd, X_pred, ncomp = 3)

# Monitor new data
X_new <- cbind(rnorm(10), rnorm(10))
Y_new <- matrix(rnorm(10 * m) + 2, 10, m)
fd_new <- fdata(Y_new, argvals = argvals)
mon <- frcc.monitor(chart, fd_new, X_new)
mon
#> SPM Monitoring Result (Phase II)
#>   Observations: 10 
#>   T2 alarms: 6 of 10 (60%) 
#>   SPE alarms: 10 of 10 (100%) 
#>   T2 UCL: 7.815 
#>   SPE UCL: 1.075 
# }
```
