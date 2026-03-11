# LOO-CV and PRESS

Computes leave-one-out cross-validation statistics.

## Usage

``` r
fregre.loo(model, data, y)
```

## Arguments

- model:

  A fitted `fregre.lm` model.

- data:

  An fdata object (the training data).

- y:

  Response vector.

## Value

A list with press, loo_r_squared, loo_residuals, leverage, and tss.
