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

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
result <- fregre.loo(fit, fd, y)
# }
```
