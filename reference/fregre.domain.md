# Domain Selection

Identifies the most important domain intervals for prediction using a
sliding window on pointwise importance.

## Usage

``` r
fregre.domain(model, window.width = 5, threshold = 0.1)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- window.width:

  Width of the sliding window (default 5).

- threshold:

  Minimum importance threshold (default 0.1).

## Value

A list with pointwise_importance, intervals, window_width, and
threshold.
