# Monitor New Functional Data (Phase II)

Projects new functional observations through a trained SPM chart and
checks whether T-squared and SPE statistics exceed control limits.

## Usage

``` r
spm.monitor(chart, newdata)
```

## Arguments

- chart:

  An object of class `spm.chart` from
  [`spm.phase1`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md).

- newdata:

  An object of class `fdata` with new observations to monitor.

## Value

An object of class `spm.monitor` with components:

- t2:

  T-squared values for new observations

- spe:

  SPE values for new observations

- t2.alarm:

  Logical vector: TRUE where T-squared exceeds UCL

- spe.alarm:

  Logical vector: TRUE where SPE exceeds UCL

- scores:

  FPC score matrix for new observations

- t2.ucl:

  T-squared control limit

- spe.ucl:

  SPE control limit

## See also

[`spm.phase1`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md)
for building the chart

## Examples

``` r
# \donttest{
# Build chart from in-control data
set.seed(1)
n <- 50; m <- 30
argvals <- seq(0, 1, length.out = m)
X <- matrix(rnorm(n * m), n, m)
fd <- fdata(X, argvals = argvals)
chart <- spm.phase1(fd, ncomp = 3)

# Monitor new data (with a shift)
X_new <- matrix(rnorm(10 * m) + 2, 10, m)
fd_new <- fdata(X_new, argvals = argvals)
mon <- spm.monitor(chart, fd_new)
mon
#> SPM Monitoring Result (Phase II)
#>   Observations: 10 
#>   T2 alarms: 0 of 10 (0%) 
#>   SPE alarms: 10 of 10 (100%) 
#>   T2 UCL: 7.815 
#>   SPE UCL: 1.862 
# }
```
