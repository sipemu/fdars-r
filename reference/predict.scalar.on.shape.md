# Predict from a Scalar-on-Shape Model

Predict from a Scalar-on-Shape Model

## Usage

``` r
# S3 method for class 'scalar.on.shape'
predict(object, newdata, ...)
```

## Arguments

- object:

  A fitted `scalar.on.shape` model.

- newdata:

  An object of class 'fdata'.

- ...:

  Ignored.

## Value

Predicted response values (numeric vector).

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- scalar.on.shape(fd, y, nbasis = 5)
pred <- predict(fit, fd[1:10, ])
# }
```
