# Accumulated Local Effects

Computes ALE for a specific FPC component.

## Usage

``` r
fregre.ale(model, data, component, n.bins = 20)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- data:

  An fdata object.

- component:

  Which FPC component (1-based).

- n.bins:

  Number of ALE bins (default 20).

## Value

A list with bin_midpoints, ale_values, bin_edges, bin_counts, and
component.
