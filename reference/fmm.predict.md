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
