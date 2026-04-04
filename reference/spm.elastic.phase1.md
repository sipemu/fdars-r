# Elastic SPM Phase I Chart

Builds an SPM control chart that is invariant to phase (time-warping)
variation. First aligns the functional data to a Karcher mean using
elastic registration, then builds separate amplitude and (optionally)
phase charts.

## Usage

``` r
spm.elastic.phase1(
  fdataobj,
  ncomp = 5,
  alpha = 0.05,
  tuning.fraction = 0.5,
  seed = 42,
  align.lambda = 0,
  monitor.phase = TRUE,
  warp.ncomp = 3
)
```

## Arguments

- fdataobj:

  An object of class `fdata` with in-control data.

- ncomp:

  Number of principal components for amplitude chart (default 5).

- alpha:

  Significance level (default 0.05).

- tuning.fraction:

  Fraction for FPCA tuning (default 0.5).

- seed:

  Random seed (default 42).

- align.lambda:

  Elastic alignment regularization (default 0.0).

- monitor.phase:

  Logical; if `TRUE`, also builds a phase chart for warping function
  monitoring (default `TRUE`).

- warp.ncomp:

  Number of components for the phase (warping) chart (default 3).

## Value

An object of class `spm.elastic.chart` with components:

- karcher.mean:

  The Karcher mean function

- mean.alignment.residual:

  Residual from Karcher mean computation

- amplitude:

  List with amplitude chart fields: eigenvalues, t2.phase1, spe.phase1,
  t2.ucl, spe.ucl, ncomp

- phase:

  List with phase chart fields (if `monitor.phase = TRUE`), or NULL

- fdataobj:

  Original fdata object

- .rust:

  Internal fields for Phase II monitoring

## See also

[`spm.elastic.monitor`](https://sipemu.github.io/fdars-r/reference/spm.elastic.monitor.md)
for Phase II elastic monitoring,
[`spm.phase1`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md)
for non-elastic Phase I

## Examples

``` r
# \donttest{
set.seed(1)
n <- 40; m <- 30
argvals <- seq(0, 1, length.out = m)
X <- matrix(rnorm(n * m), n, m)
fd <- fdata(X, argvals = argvals)

chart <- spm.elastic.phase1(fd, ncomp = 3, warp.ncomp = 2)
chart
#> Elastic SPM Control Chart (Phase I)
#>   Components (amplitude): 3 
#>   Alpha: 0.05 
#>   Amplitude T2 UCL: 7.815 
#>   Amplitude SPE UCL: 0.7933 
#>   Phase monitoring: enabled
#>   Phase components: 2 
#>   Phase T2 UCL: 5.991 
#>   Phase SPE UCL: 0.003706 
#>   Observations: 40 
#>   Grid points: 30 
# }
```
