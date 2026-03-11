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
