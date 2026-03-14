# Build Functional Regression Control Chart (Phase I)

Monitors a functional response after adjusting for known scalar
covariates using function-on-scalar regression (FOSR). The residuals are
then monitored via FPCA-based T-squared and SPE statistics.

## Usage

``` r
frcc.phase1(
  fdataobj,
  predictors,
  ncomp = 5,
  fosr.lambda = 1e-04,
  alpha = 0.05,
  tuning.fraction = 0.5,
  seed = 42
)
```

## Arguments

- fdataobj:

  An object of class `fdata` (functional response).

- predictors:

  A matrix of scalar predictors (n x p).

- ncomp:

  Number of principal components for residual FPCA (default 5).

- fosr.lambda:

  FOSR smoothing parameter (default 1e-4).

- alpha:

  Significance level (default 0.05).

- tuning.fraction:

  Fraction of data for tuning (default 0.5).

- seed:

  Random seed (default 42).

## Value

An object of class `frcc.chart` with components:

- eigenvalues:

  Eigenvalues from residual FPCA

- t2.ucl:

  T-squared control limit

- spe.ucl:

  SPE control limit

- ncomp:

  Number of components used

- fdataobj:

  Original fdata object

- predictors:

  Original predictor matrix

- .rust:

  Internal fields for Phase II monitoring

## See also

[`frcc.monitor`](https://sipemu.github.io/fdars-r/reference/frcc.monitor.md)
for Phase II monitoring,
[`spm.phase1`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md)
for monitoring without covariates

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
chart
#> Functional Regression Control Chart (Phase I)
#>   Components: 3 
#>   Alpha: 0.05 
#>   T2 UCL: 7.815 
#>   SPE UCL: 1.171 
#>   Observations: 60 
#>   Predictors: 2 
# }
```
