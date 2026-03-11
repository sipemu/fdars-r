# Conditional Permutation Importance

Permutation importance respecting correlations between components.

## Usage

``` r
fregre.conditional.importance(
  model,
  data,
  y,
  n.bins = 5,
  n.perm = 100,
  seed = NULL
)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- data:

  An fdata object.

- y:

  Response vector.

- n.bins:

  Number of conditional bins (default 5).

- n.perm:

  Number of permutations (default 100).

- seed:

  Random seed.

## Value

A list with importance, baseline_metric, and unconditional_importance.
