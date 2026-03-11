# Influence Diagnostics

Computes leverage, Cook's distance, and related diagnostics.

## Usage

``` r
fregre.influence(model, data)
```

## Arguments

- model:

  A fitted `fregre.lm` model.

- data:

  An fdata object (the training data).

## Value

A list with leverage, cooks_distance, p, and mse.
