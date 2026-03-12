# Prediction Intervals

Computes prediction intervals for new observations.

## Usage

``` r
fregre.prediction.interval(model, train.data, new.data, confidence = 0.95)
```

## Arguments

- model:

  A fitted `fregre.lm` model.

- train.data:

  An fdata object (training data).

- new.data:

  An fdata object (new observations for prediction).

- confidence:

  Confidence level (default 0.95).

## Value

A list with predictions, lower, upper, prediction_se, confidence_level,
t_critical, and residual_se.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
new_fd <- fdata(matrix(rnorm(50), nrow = 5), argvals = seq(0, 1, length.out = 10))
result <- fregre.prediction.interval(fit, fd, new_fd)
# }
```
