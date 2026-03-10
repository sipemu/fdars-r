# Create Stratified or Random Fold Assignments

Create Stratified or Random Fold Assignments

## Usage

``` r
.create_folds(y, kfold, type, stratified, seed)
```

## Arguments

- y:

  Response vector.

- kfold:

  Number of folds.

- type:

  `"regression"` or `"classification"`.

- stratified:

  Logical.

- seed:

  Optional seed.

## Value

Integer vector of fold assignments.
