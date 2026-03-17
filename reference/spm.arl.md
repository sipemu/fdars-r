# Estimate Average Run Length (ARL)

Estimates the Average Run Length (ARL) for a T-squared or SPE control
chart via Monte Carlo simulation. In-control ARL
(`shift.magnitude = NULL`) measures the expected number of observations
before a false alarm. Out-of-control ARL measures detection speed under
a mean shift.

## Usage

``` r
spm.arl(
  chart,
  type = c("t2", "spe"),
  shift.magnitude = NULL,
  n.sim = 10000,
  max.rl = 5000,
  seed = 42
)
```

## Arguments

- chart:

  An object of class `spm.chart` from
  [`spm.phase1`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md).

- type:

  Which statistic: `"t2"` or `"spe"`.

- shift.magnitude:

  For out-of-control ARL: numeric vector of per-component shifts (length
  = ncomp). If `NULL`, computes in-control ARL.

- n.sim:

  Number of Monte Carlo simulations (default 10000).

- max.rl:

  Maximum run length per simulation (default 5000).

- seed:

  Random seed (default 42).

## Value

An object of class `spm.arl` with components:

- arl:

  Estimated ARL

- std.dev:

  Standard deviation of run lengths

- median.rl:

  Median run length

- run.lengths:

  Integer vector of all simulated run lengths

- type:

  Which statistic was used

- in.control:

  Logical: TRUE if in-control ARL

## See also

[`spm.arl.ewma`](https://sipemu.github.io/fdars-r/reference/spm.arl.ewma.md)
for EWMA-based ARL

## Examples

``` r
# \donttest{
set.seed(1)
n <- 50; m <- 30
argvals <- seq(0, 1, length.out = m)
X <- matrix(rnorm(n * m), n, m)
fd <- fdata(X, argvals = argvals)
chart <- spm.phase1(fd, ncomp = 3)

# In-control ARL
arl0 <- spm.arl(chart, type = "t2", n.sim = 1000)
arl0
#> SPM Average Run Length (In-Control)
#>   Statistic: T2 
#>   ARL: 19.62 
#>   Std. Dev.: 19.02 
#>   Median RL: 14 
#>   Simulations: 1000 
# }
```
