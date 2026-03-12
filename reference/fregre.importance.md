# FPC Permutation Importance

Measures importance of each FPC by permuting scores and measuring the
increase in prediction error.

## Usage

``` r
fregre.importance(model, data, y, n.perm = 100, seed = NULL)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- data:

  An fdata object (the training data).

- y:

  Response vector.

- n.perm:

  Number of permutations (default 100).

- seed:

  Random seed.

## Value

A list with importance, baseline_metric, and permuted_metric.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
result <- fregre.importance(fit, fd, y)
# }
```
