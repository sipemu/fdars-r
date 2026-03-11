# Variance Inflation Factors

Computes VIF for FPC scores to detect multicollinearity.

## Usage

``` r
fregre.vif(model, data)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- data:

  An fdata object (the training data).

## Value

A list with vif, labels, mean_vif, n_moderate, and n_severe.
