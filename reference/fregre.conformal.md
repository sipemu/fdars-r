# Conformal Prediction

Computes distribution-free prediction intervals using conformal
inference.

## Usage

``` r
fregre.conformal(
  model,
  train.data,
  train.y,
  test.data,
  cal.fraction = 0.2,
  alpha = 0.1,
  seed = NULL
)
```

## Arguments

- model:

  A fitted `fregre.lm` model.

- train.data:

  An fdata object (training data).

- train.y:

  Training response vector.

- test.data:

  An fdata object (test data).

- cal.fraction:

  Fraction of training data for calibration (default 0.2).

- alpha:

  Miscoverage level (default 0.1 for 90% intervals).

- seed:

  Random seed.

## Value

A list with predictions, lower, upper, residual_quantile, and coverage.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
new_fd <- fdata(matrix(rnorm(50), nrow = 5), argvals = seq(0, 1, length.out = 10))
result <- fregre.conformal(fit, fd, y, new_fd)
# }
```
