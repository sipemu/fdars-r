# Cross-Validated Functional Classification

Evaluates classification error rate using k-fold cross-validation.

## Usage

``` r
fclassif.cv(
  fdataobj,
  y,
  method = "lda",
  covariates = NULL,
  ncomp = 3,
  nfold = 10,
  seed = NULL
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- y:

  Integer vector of class labels.

- method:

  Classification method (default "lda").

- covariates:

  Optional scalar covariates matrix.

- ncomp:

  Number of FPC components (default 3).

- nfold:

  Number of CV folds (default 10).

- seed:

  Random seed for fold assignment.

## Value

An object of class 'fclassif.cv' with components:

- error.rate:

  Mean error rate across folds

- fold.errors:

  Per-fold error rates

- best.ncomp:

  Best ncomp if tuned
