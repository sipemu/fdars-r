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

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
result <- fregre.ale(fit, fd, component = 1)
# }
```
