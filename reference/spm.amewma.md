# Adaptive MEWMA (AMEWMA) Monitoring for Functional Data

Applies Adaptive Multivariate EWMA monitoring where the smoothing
parameter adapts over time based on the magnitude of recent deviations.
This provides good detection performance for both small and large
shifts.

## Usage

``` r
spm.amewma(
  chart,
  newdata,
  lambda.min = 0.05,
  lambda.max = 0.95,
  lambda.init = 0.2,
  eta = 0.1,
  ncomp = NULL,
  alpha = 0.05
)
```

## Arguments

- chart:

  An object of class `spm.chart` from
  [`spm.phase1`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md).

- newdata:

  An object of class `fdata` with sequential observations.

- lambda.min:

  Minimum smoothing parameter (default 0.05).

- lambda.max:

  Maximum smoothing parameter (default 0.95).

- lambda.init:

  Initial smoothing parameter (default 0.2).

- eta:

  Learning rate for lambda adaptation (default 0.1).

- ncomp:

  Number of components (default: same as chart).

- alpha:

  Significance level (default 0.05).

## Value

An object of class `spm.amewma` with components:

- smoothed.scores:

  Matrix of adaptively smoothed score vectors

- t2.statistic:

  Numeric vector of T-squared statistics on smoothed scores

- lambda.t:

  Numeric vector of adaptive lambda values over time

- ucl:

  Control limit

- alarm:

  Logical vector: TRUE where statistic exceeds UCL

## See also

[`spm.mewma`](https://sipemu.github.io/fdars-r/reference/spm.mewma.md)
for fixed-lambda MEWMA,
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

X_new <- matrix(rnorm(30 * m) + 0.5, 30, m)
fd_new <- fdata(X_new, argvals = argvals)
amewma <- spm.amewma(chart, fd_new)
amewma
#> SPM Adaptive MEWMA (AMEWMA) Monitoring
#>   Observations: 30 
#>   Lambda range: [0.05, 0.95]
#>   Lambda init: 0.2 
#>   Eta: 0.1 
#>   UCL: 7.815 
#>   Alarms: 1 of 30 (3.33%) 
# }
```
