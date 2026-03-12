# FPC-based Function-on-Scalar Regression

Fits a function-on-scalar regression using functional principal
components.

## Usage

``` r
fosr.fpc(fdataobj, predictors, ncomp = 3)
```

## Arguments

- fdataobj:

  An object of class 'fdata' (functional response).

- predictors:

  A matrix of scalar predictors (n x p).

- ncomp:

  Number of FPC components (default 3).

## Value

An object of class 'fosr'.

## Examples

``` r
# \donttest{
Y <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
X <- cbind(rnorm(50), rnorm(50))
fit <- fosr.fpc(Y, predictors = X, ncomp = 2)
fit
#> Function-on-Scalar Regression
#> =============================
#>   Number of observations: 50 
#>   Number of predictors: 2 
#>   Evaluation points: 10 
#>   R-squared: 0.0037 
#>   FPC components: 2 
# }
```
