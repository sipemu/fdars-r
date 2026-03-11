# Functional Partial Dependence Plot

Computes the partial dependence of the response on a specific FPC score,
along with ICE (Individual Conditional Expectation) curves.

## Usage

``` r
fregre.pdp(model, data, component, n.grid = 20)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- data:

  An fdata object (the training data).

- component:

  Which FPC component (1-based).

- n.grid:

  Number of grid points for the PDP (default 20).

## Value

A list with grid_values, pdp_curve, ice_curves, and component.
