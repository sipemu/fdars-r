# Functional Linear Model (FPC-based)

Fits a functional linear model using FPC regression with the Rust
backend. Supports optional scalar covariates.

## Usage

``` r
fregre.lm(fdataobj, y, scalar.covariates = NULL, ncomp = NULL)
```

## Arguments

- fdataobj:

  An object of class 'fdata' (functional covariate).

- y:

  Response vector (scalar).

- scalar.covariates:

  Optional matrix of scalar covariates (n x p).

- ncomp:

  Number of FPC components to use. If NULL, uses min(n-1, 15).

## Value

A fitted regression object of class 'fregre.lm' with components:

- intercept:

  Intercept term

- beta.t:

  Functional coefficient beta(t) as fdata

- gamma:

  Scalar covariate coefficients

- fitted.values:

  Fitted values

- residuals:

  Residuals

- r.squared:

  R-squared

- r.squared.adj:

  Adjusted R-squared

- ncomp:

  Number of FPC components used

- gcv:

  GCV criterion value

## See also

[`fregre.pc`](https://sipemu.github.io/fdars-r/reference/fregre.pc.md)
for the R-native FPC regression,
[`fregre.lm.cv`](https://sipemu.github.io/fdars-r/reference/fregre.lm.cv.md)
for cross-validated component selection

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
result <- fregre.lm(fd, y, ncomp = 3)
result$r.squared
#> [1] 0.1393605
# }
```
