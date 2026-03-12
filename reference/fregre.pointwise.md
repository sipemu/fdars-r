# Pointwise Importance

Computes importance of each domain point for prediction.

## Usage

``` r
fregre.pointwise(model)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

## Value

A list with importance, importance_normalized, component_importance, and
score_variance.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
result <- fregre.pointwise(fit)
# }
```
