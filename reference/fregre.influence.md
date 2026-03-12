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

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
result <- fregre.influence(fit, fd)
# }
```
