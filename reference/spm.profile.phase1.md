# Profile Monitoring Phase I

Builds a profile monitoring chart for functional data with scalar
covariates. Uses function-on-scalar regression (FOSR) to model the
relationship, then applies FPCA to the regression coefficients across
sliding windows.

## Usage

``` r
spm.profile.phase1(
  fdataobj,
  predictors,
  ncomp = 3,
  alpha = 0.05,
  fosr.lambda = 1e-04,
  window.size = 20,
  step.size = 1
)
```

## Arguments

- fdataobj:

  An object of class `fdata` (functional response).

- predictors:

  A matrix of scalar predictors (n x p).

- ncomp:

  Number of principal components for beta FPCA (default 3).

- alpha:

  Significance level (default 0.05).

- fosr.lambda:

  FOSR smoothing parameter (default 1e-4).

- window.size:

  Window size for sliding-window FOSR (default 20).

- step.size:

  Step size between windows (default 1).

## Value

An object of class `spm.profile.chart` with components:

- eigenvalues:

  Eigenvalues from beta FPCA

- t2.ucl:

  T-squared control limit

- t2.description:

  Description of the limit

- lag1.autocorrelation:

  Lag-1 autocorrelation of T-squared values

- effective.n.windows:

  Effective number of independent windows

- fdataobj:

  Original fdata object

- predictors:

  Original predictor matrix

- .rust:

  Internal fields for Phase II monitoring

## See also

[`spm.profile.monitor`](https://sipemu.github.io/fdars-r/reference/spm.profile.monitor.md)
for Phase II monitoring

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
chart
#> SPM Profile Monitoring Chart (Phase I)
#>   Components: 2 
#>   Alpha: 0.05 
#>   T2 UCL: 5.991 
#>   Observations: 80 
#>   Predictors: 2 
#>   Lag-1 autocorrelation: 0.8858 
#>   Effective windows: 24 
# }
```
