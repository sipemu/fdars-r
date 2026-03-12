# Function-on-Scalar Regression

Functions for regression where the response is a function and predictors
are scalar variables. Function-on-Scalar Regression (Penalized)

## Usage

``` r
fosr(fdataobj, predictors, lambda = 0)
```

## Arguments

- fdataobj:

  An object of class 'fdata' (functional response).

- predictors:

  A matrix of scalar predictors (n x p).

- lambda:

  Smoothing/penalty parameter (default 0).

## Value

An object of class 'fosr' with components:

- intercept:

  Intercept function mu(t) as fdata

- beta:

  Coefficient functions beta_j(t) as fdata

- fitted:

  Fitted functional values as fdata

- residuals:

  Residual functions as fdata

- r.squared:

  Global R-squared

- r.squared.t:

  Pointwise R-squared

- gcv:

  GCV criterion value

- lambda:

  Smoothing parameter used

## Details

Fits a penalized function-on-scalar regression model where the response
is a function Y_i(t) and predictors are scalar:
`Y_i(t) = mu(t) + sum_j x_ij beta_j(t) + eps_i(t)`.

## See also

[`fosr.fpc`](https://sipemu.github.io/fdars-r/reference/fosr.fpc.md) for
FPC-based FOSR,
[`fanova`](https://sipemu.github.io/fdars-r/reference/fanova.md) for
functional ANOVA

## Examples

``` r
# \donttest{
# Functional response: 50 curves observed at 10 time points
Y <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
# Two scalar predictors
X <- cbind(rnorm(50), rnorm(50))
fit <- fosr(Y, predictors = X, lambda = 1)
fit
#> Function-on-Scalar Regression
#> =============================
#>   Number of observations: 50 
#>   Number of predictors: 2 
#>   Evaluation points: 10 
#>   R-squared: 0.0302 
#>   Lambda: 1 
# }
```
