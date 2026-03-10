# Unified K-Fold Cross-Validation for Functional Data

Performs k-fold cross-validation for any functional data model,
producing out-of-fold (OOF) predictions where each observation is
predicted exactly once (when it is in the test fold). Supports both
regression and classification, with optional stratified fold assignment.
Repeated cross-validation (`nrep > 1`) runs multiple random fold
partitions to assess prediction variability.

## Usage

``` r
cv.fdata(
  fdataobj,
  y,
  fit.fn,
  predict.fn = NULL,
  kfold = 10,
  nrep = 1,
  type = c("regression", "classification"),
  stratified = TRUE,
  seed = NULL,
  ...
)
```

## Arguments

- fdataobj:

  An `fdata` object containing functional predictors.

- y:

  Response vector: numeric for regression, factor (or coerced to factor)
  for classification.

- fit.fn:

  A function with signature `function(fdataobj, y, ...)` that returns a
  fitted model object.

- predict.fn:

  Optional prediction function with signature
  `function(model, newdata, ...)`. If `NULL` (default), uses
  `predict(model, newdata)`.

- kfold:

  Number of folds (default 10).

- nrep:

  Number of repetitions (default 1). When `nrep > 1`, the entire k-fold
  procedure is repeated with different random fold assignments,
  producing an `n x nrep` matrix of predictions and per-observation
  variability estimates.

- type:

  One of `"regression"` or `"classification"`. Auto-detected from `y` if
  it is a factor.

- stratified:

  Logical; whether to stratify fold assignments so each fold has a
  similar distribution of `y` (default `TRUE`).

- seed:

  Optional integer seed for reproducibility of fold assignments.

- ...:

  Additional arguments passed through to `fit.fn`.

## Value

An object of class `"cv.fdata"` with components:

- oof.predictions:

  Out-of-fold predictions (numeric vector for regression, factor for
  classification). When `nrep > 1`, aggregated across repetitions (mean
  for regression, majority vote for classification).

- oof.probabilities:

  For classification only: matrix of class probabilities if `predict.fn`
  returns them; otherwise `NULL`.

- folds:

  Integer vector of fold assignments (from rep 1 for backward
  compatibility).

- fold.models:

  List of per-fold fitted model objects. When `nrep > 1`, a list of
  `nrep` lists, each of length `kfold`.

- metrics:

  Named list of overall performance metrics (computed on aggregated
  predictions when `nrep > 1`).

- fold.metrics:

  Data frame of per-fold metrics (from rep 1 for backward
  compatibility).

- call:

  The matched call.

- type:

  Character: `"regression"` or `"classification"`.

- kfold:

  Number of folds used.

- y:

  The response vector.

- nrep:

  Number of repetitions (only when `nrep > 1`).

- oof.matrix:

  Matrix (`n x nrep`) of per-repetition predictions (only when
  `nrep > 1`).

- oof.sd:

  Numeric vector of per-observation prediction SD across repetitions
  (only when `nrep > 1`). For classification, the proportion of reps
  disagreeing with the majority vote.

- folds.matrix:

  Integer matrix (`n x nrep`) of fold assignments per repetition (only
  when `nrep > 1`).

- rep.metrics:

  Data frame with one row per repetition containing per-rep metrics
  (only when `nrep > 1`).

- metrics.summary:

  List with mean and sd of each metric across repetitions (only when
  `nrep > 1`).

## Details

The `fit.fn` can internally perform nested cross-validation for
hyperparameter tuning. For example, passing a wrapper around
[`fregre.pc.cv`](https://sipemu.github.io/fdars-r/reference/fregre.pc.cv.md)
as `fit.fn` gives proper nested CV: outer folds produce unbiased OOF
predictions while inner CV selects optimal parameters on training data
only.

Each fold is wrapped in `tryCatch`: if a fold fails, its predictions are
set to `NA` and a warning is issued, but the remaining folds continue.

When `nrep > 1`, each repetition uses a different random fold partition.
If `seed` is provided, deterministic per-rep seeds are derived from it
for full reproducibility.

## Examples

``` r
# Simple regression CV with fixed hyperparameters
set.seed(1)
t <- seq(0, 1, length.out = 30)
X <- matrix(0, 40, 30)
for (i in 1:40) X[i, ] <- sin(2*pi*t * i/40) + rnorm(30, sd = 0.1)
y <- rowMeans(X) + rnorm(40, sd = 0.1)
fd <- fdata(X, argvals = t)

result <- cv.fdata(fd, y,
  fit.fn = function(fd, y, ...) fregre.pc(fd, y, ncomp = 3),
  kfold = 5, seed = 42)
print(result)
#> K-Fold Cross-Validation (cv.fdata)
#>   Type: regression 
#>   Folds: 5 
#>   Observations: 40  
#> 
#> Overall metrics:
#>   RMSE: 0.1117 
#>   MAE:  0.08839 
#>   R2:   0.8551 
#> 
#> Per-fold RMSE range: [0.05177, 0.1356]
```
