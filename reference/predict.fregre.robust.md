# Predict from Robust Functional Regression

Predict from Robust Functional Regression

## Usage

``` r
# S3 method for class 'fregre.robust'
predict(object, newdata, scalar.covariates = NULL, ...)
```

## Arguments

- object:

  A fitted 'fregre.robust' model.

- newdata:

  An object of class 'fdata' (new curves).

- scalar.covariates:

  Optional new scalar covariates.

- ...:

  Additional arguments (ignored).

## Value

Predicted response values.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.l1(fd, y, ncomp = 3)
pred <- predict(fit, fd[1:5, ])
# }
```
