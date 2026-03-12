# Functional Saliency Map

Computes a saliency map showing which domain regions most influence each
observation's prediction.

## Usage

``` r
fregre.saliency(model, data = NULL)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- data:

  An fdata object (the training data, used for lm models).

## Value

A list with saliency_map (matrix) and mean_absolute_saliency.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
result <- fregre.saliency(fit, fd)
# }
```
