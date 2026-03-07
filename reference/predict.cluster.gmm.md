# Predict Cluster Membership for New Functional Data

Assigns new functional observations to clusters using a fitted GMM.

## Usage

``` r
# S3 method for class 'cluster.gmm'
predict(object, newdata, new.covariates = NULL, ...)
```

## Arguments

- object:

  A fitted 'cluster.gmm' object.

- newdata:

  An fdata object of new curves.

- new.covariates:

  Optional matrix of covariates for new data.

- ...:

  Additional arguments (ignored).

## Value

A list with `cluster` (assignments) and `membership` (probabilities).
