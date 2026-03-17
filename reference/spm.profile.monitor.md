# Monitor New Data Against Profile Chart (Phase II)

Projects new functional data and scalar predictors through a trained
profile monitoring chart. Computes sliding-window FOSR betas, projects
onto the reference beta FPCA, and evaluates T-squared statistics.

## Usage

``` r
spm.profile.monitor(chart, newdata, new.predictors)
```

## Arguments

- chart:

  An object of class `spm.profile.chart` from
  [`spm.profile.phase1`](https://sipemu.github.io/fdars-r/reference/spm.profile.phase1.md).

- newdata:

  An object of class `fdata` with new functional response.

- new.predictors:

  A matrix of new scalar predictors.

## Value

An object of class `spm.monitor` with components:

- betas:

  Matrix of estimated beta coefficients per window

- t2:

  T-squared values for each window

- t2.alarm:

  Logical: TRUE where T-squared exceeds UCL

- beta.scores:

  Beta FPC scores for each window

- t2.ucl:

  T-squared control limit

## See also

[`spm.profile.phase1`](https://sipemu.github.io/fdars-r/reference/spm.profile.phase1.md)
for building the chart

## Examples

``` r
# \donttest{
set.seed(1)
n <- 80; m <- 20
argvals <- seq(0, 1, length.out = m)
X_pred <- cbind(rnorm(n), rnorm(n))
Y <- matrix(rnorm(n * m), n, m)
fd <- fdata(Y, argvals = argvals)
chart <- spm.profile.phase1(fd, X_pred, ncomp = 2, window.size = 15)

# Monitor new data
X_new <- cbind(rnorm(30), rnorm(30))
Y_new <- matrix(rnorm(30 * m) + 1, 30, m)
fd_new <- fdata(Y_new, argvals = argvals)
mon <- spm.profile.monitor(chart, fd_new, X_new)
mon
#> SPM Monitoring Result (Phase II)
#>   Observations: 16 
#>   T2 alarms: 0 of 16 (0%) 
#>   SPE alarms: 0 of 16 (0%) 
#>   T2 UCL: 5.991 
#>   SPE UCL: NULL 
# }
```
