# Build Univariate SPM Control Chart (Phase I)

Builds an FPCA-based control chart from in-control functional data. The
data is split into a tuning set (for FPCA) and a calibration set (for
control limits). Computes Hotelling T-squared and SPE statistics.

## Usage

``` r
spm.phase1(fdataobj, ncomp = 5, alpha = 0.05, tuning.fraction = 0.5, seed = 42)
```

## Arguments

- fdataobj:

  An object of class `fdata` containing in-control functional data (at
  least 4 observations).

- ncomp:

  Number of principal components (default 5).

- alpha:

  Significance level for control limits (default 0.05).

- tuning.fraction:

  Fraction of data used for FPCA tuning (default 0.5).

- seed:

  Random seed for train/calibration split (default 42).

## Value

An object of class `spm.chart` with components:

- eigenvalues:

  Eigenvalues from FPCA

- t2.phase1:

  T-squared values for calibration set

- spe.phase1:

  SPE values for calibration set

- t2.ucl:

  T-squared upper control limit

- spe.ucl:

  SPE upper control limit

- ncomp:

  Number of components used

- t2.description:

  Description of T-squared limit

- spe.description:

  Description of SPE limit

- fdataobj:

  Original fdata object

- .rust:

  Internal fields for Phase II monitoring

## See also

[`spm.monitor`](https://sipemu.github.io/fdars-r/reference/spm.monitor.md)
for Phase II monitoring,
[`spm.ewma`](https://sipemu.github.io/fdars-r/reference/spm.ewma.md) for
EWMA-based monitoring,
[`frcc.phase1`](https://sipemu.github.io/fdars-r/reference/frcc.phase1.md)
for regression-adjusted monitoring

## Examples

``` r
# \donttest{
# Simulate in-control functional data
set.seed(1)
n <- 50
m <- 30
argvals <- seq(0, 1, length.out = m)
X <- matrix(rnorm(n * m), n, m)
fd <- fdata(X, argvals = argvals)

# Build Phase I chart
chart <- spm.phase1(fd, ncomp = 3)
chart
#> SPM Control Chart (Phase I)
#>   Components: 3 
#>   Alpha: 0.05 
#>   T2 UCL: 7.815 
#>   SPE UCL: 1.811 
#>   Observations: 50 
#>   Grid points: 30 
#>   Eigenvalues: 0.158, 0.122, 0.107 
# }
```
