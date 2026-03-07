# Cross-Validation for FPC Component Selection (fregre.lm)

Uses the Rust backend to select the optimal number of FPC components via
k-fold cross-validation for the functional linear model.

## Usage

``` r
fregre.lm.cv(fdataobj, y, scalar.covariates = NULL, k.range = NULL, nfold = 10)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- y:

  Response vector.

- scalar.covariates:

  Optional matrix of scalar covariates.

- k.range:

  Range of FPC component counts to try.

- nfold:

  Number of CV folds (default 10).

## Value

A list with `optimal.k`, `cv.errors`, and `model`.
