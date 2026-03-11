# Regression Depth

Computes regression depth for the fitted model and individual
observation score depths.

## Usage

``` r
fregre.depth(model, data, y, n.boot = 100, seed = NULL)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- data:

  An fdata object.

- y:

  Response vector.

- n.boot:

  Number of bootstrap resamples (default 100).

- seed:

  Random seed.

## Value

A list with beta_depth, score_depths, mean_score_depth, and
n_boot_success.
