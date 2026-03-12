# Predict from Functional Mixed Model

Predict population-level (fixed effects only) functional responses for
new covariate values.

## Usage

``` r
fmm.predict(object, new.covariates = NULL, ...)
```

## Arguments

- object:

  A fitted 'fmm' object.

- new.covariates:

  Matrix of new covariate values. If NULL, returns the overall mean
  function.

- ...:

  Additional arguments (ignored).

## Value

An fdata object of predicted functional values.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
subject <- rep(1:10, each = 5)
x <- cbind(rnorm(50))
fit <- fmm(fd, subject.ids = subject, covariates = x)
pred <- fmm.predict(fit, new.covariates = cbind(c(0, 1)))
# }
```
