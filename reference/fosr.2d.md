# 2D Function-on-Scalar Regression

Fits a penalized regression model where the response is a 2D surface
Y_i(s,t) and predictors are scalar.

## Usage

``` r
fosr.2d(fdataobj, predictors, argvals.s, argvals.t, lambda.s = 0, lambda.t = 0)
```

## Arguments

- fdataobj:

  An fdata object where each row is a flattened 2D surface (m1 \* m2
  grid points).

- predictors:

  A matrix of scalar predictors (n x p).

- argvals.s:

  Grid points along the s dimension.

- argvals.t:

  Grid points along the t dimension.

- lambda.s:

  Smoothing penalty in s direction (default 0).

- lambda.t:

  Smoothing penalty in t direction (default 0).

## Value

An object of class 'fosr.2d' with components:

- intercept:

  Intercept surface as fdata

- beta:

  Coefficient surfaces as fdata

- fitted:

  Fitted surfaces as fdata

- residuals:

  Residual surfaces as fdata

- r.squared:

  Global R-squared

- r.squared.pointwise:

  Pointwise R-squared values

- lambda.s:

  Penalty in s direction

- lambda.t:

  Penalty in t direction

- gcv:

  GCV criterion

- grid:

  List with argvals.s, argvals.t, m1, m2

## Examples

``` r
# \donttest{
n <- 30; m1 <- 5; m2 <- 5
fd <- fdata(matrix(rnorm(n * m1 * m2), n, m1 * m2))
X <- matrix(rnorm(n * 2), n, 2)
fit <- fosr.2d(fd, X, seq(0, 1, length.out = m1), seq(0, 1, length.out = m2))
# }
```
