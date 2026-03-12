# Variance Inflation Factors

Computes VIF for FPC scores to detect multicollinearity.

## Usage

``` r
fregre.vif(model, data)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- data:

  An fdata object (the training data).

## Value

A list with vif, labels, mean_vif, n_moderate, and n_severe.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
result <- fregre.vif(fit, fd)
# }
```
