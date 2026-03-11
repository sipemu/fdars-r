# Friedman H-Statistic

Measures the interaction strength between two FPC components using the
Friedman H-statistic.

## Usage

``` r
fregre.friedman(model, data, component.j, component.k, n.grid = 20)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- data:

  An fdata object.

- component.j:

  First component index (1-based).

- component.k:

  Second component index (1-based).

- n.grid:

  Grid size for PDP computation (default 20).

## Value

A list with component_j, component_k, h_squared, grid_j, grid_k, and
pdp_2d.
