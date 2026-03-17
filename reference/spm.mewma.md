# MEWMA Monitoring for Functional Data

Applies Multivariate Exponentially Weighted Moving Average (MEWMA)
monitoring to sequential functional data. MEWMA smooths the FPC score
vectors over time, then computes a Hotelling-type statistic on the
smoothed scores.

## Usage

``` r
spm.mewma(chart, newdata, lambda = 0.2, ncomp = NULL, alpha = 0.05)
```

## Arguments

- chart:

  An object of class `spm.chart` from
  [`spm.phase1`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md).

- newdata:

  An object of class `fdata` with sequential observations.

- lambda:

  MEWMA smoothing parameter in (0, 1\] (default 0.2).

- ncomp:

  Number of components (default: same as chart).

- alpha:

  Significance level (default 0.05).

## Value

An object of class `spm.mewma` with components:

- smoothed.scores:

  Matrix of MEWMA-smoothed score vectors

- mewma.statistic:

  Numeric vector of MEWMA statistics

- ucl:

  Control limit

- alarm:

  Logical vector: TRUE where MEWMA exceeds UCL

- spe:

  SPE values

- spe.limit:

  SPE control limit

- spe.alarm:

  Logical: TRUE where SPE exceeds limit

## See also

[`spm.amewma`](https://sipemu.github.io/fdars-r/reference/spm.amewma.md)
for adaptive MEWMA,
[`spm.ewma`](https://sipemu.github.io/fdars-r/reference/spm.ewma.md) for
univariate EWMA

## Examples

``` r
# \donttest{
set.seed(1)
n <- 50; m <- 30
argvals <- seq(0, 1, length.out = m)
X <- matrix(rnorm(n * m), n, m)
fd <- fdata(X, argvals = argvals)
chart <- spm.phase1(fd, ncomp = 3)

X_new <- matrix(rnorm(20 * m) + 0.5, 20, m)
fd_new <- fdata(X_new, argvals = argvals)
mewma <- spm.mewma(chart, fd_new, lambda = 0.2)
mewma
#> SPM MEWMA Monitoring
#>   Observations: 20 
#>   Lambda: 0.2 
#>   UCL: 7.815 
#>   MEWMA alarms: 1 of 20 (5%) 
#>   SPE alarms: 1 of 20 (5%) 
# }
```
