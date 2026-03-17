# Estimate ARL for EWMA-T-Squared Chart

Estimates the in-control Average Run Length for an EWMA-smoothed
T-squared chart via Monte Carlo simulation.

## Usage

``` r
spm.arl.ewma(
  chart,
  lambda = 0.2,
  ucl = NULL,
  n.sim = 10000,
  max.rl = 5000,
  seed = 42
)
```

## Arguments

- chart:

  An object of class `spm.chart` from
  [`spm.phase1`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md).

- lambda:

  EWMA smoothing parameter (default 0.2).

- ucl:

  Upper control limit for the EWMA chart. If `NULL`, uses
  `chart$t2.ucl`.

- n.sim:

  Number of simulations (default 10000).

- max.rl:

  Maximum run length per simulation (default 5000).

- seed:

  Random seed (default 42).

## Value

An object of class `spm.arl` with the same structure as
[`spm.arl`](https://sipemu.github.io/fdars-r/reference/spm.arl.md).

## See also

[`spm.arl`](https://sipemu.github.io/fdars-r/reference/spm.arl.md) for
standard ARL,
[`spm.ewma`](https://sipemu.github.io/fdars-r/reference/spm.ewma.md) for
EWMA monitoring

## Examples

``` r
# \donttest{
set.seed(1)
n <- 50; m <- 30
argvals <- seq(0, 1, length.out = m)
X <- matrix(rnorm(n * m), n, m)
fd <- fdata(X, argvals = argvals)
chart <- spm.phase1(fd, ncomp = 3)

arl_ewma <- spm.arl.ewma(chart, lambda = 0.2, n.sim = 1000)
arl_ewma
#> SPM Average Run Length (In-Control)
#>   Statistic: EWMA.T2 
#>   ARL: 37.92 
#>   Std. Dev.: 34.76 
#>   Median RL: 27 
#>   Simulations: 1000 
# }
```
