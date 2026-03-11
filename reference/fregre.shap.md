# FPC SHAP Values

Computes Shapley values for the FPC scores.

## Usage

``` r
fregre.shap(model, data, n.samples = 100, seed = NULL)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- data:

  An fdata object (the training data).

- n.samples:

  Number of samples for logistic SHAP (default 100).

- seed:

  Random seed (used for logistic models).

## Value

A list with values (SHAP matrix), base_value, and mean_scores.
