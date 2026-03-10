# Run One CV Repetition

Run One CV Repetition

## Usage

``` r
.run_one_cv(fdataobj, y, fit.fn, predict.fn, folds, kfold, type, ...)
```

## Arguments

- fdataobj:

  An fdata object.

- y:

  Response vector.

- fit.fn:

  Fitting function.

- predict.fn:

  Prediction function or NULL.

- folds:

  Integer vector of fold assignments.

- kfold:

  Number of folds.

- type:

  "regression" or "classification".

- ...:

  Additional arguments passed to fit.fn.

## Value

List with oof (predictions vector) and fold.models.
