# Predict from Function-on-Scalar Regression

Predict from Function-on-Scalar Regression

## Usage

``` r
# S3 method for class 'fosr'
predict(object, new.predictors = NULL, ...)
```

## Arguments

- object:

  A fitted 'fosr' object.

- new.predictors:

  Matrix of new scalar predictors. If NULL, returns fitted.

- ...:

  Additional arguments (ignored).

## Value

An fdata object of predicted functional values.
