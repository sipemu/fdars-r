# DFBETAS and DFFITS

Computes DFBETAS (influence on each coefficient) and DFFITS (influence
on fitted values) for each observation.

## Usage

``` r
fregre.dfbetas(model, data)
```

## Arguments

- model:

  A fitted `fregre.lm` model.

- data:

  An fdata object (the training data).

## Value

A list with dfbetas (matrix), dffits, studentized_residuals,
dfbetas_cutoff, and dffits_cutoff.
