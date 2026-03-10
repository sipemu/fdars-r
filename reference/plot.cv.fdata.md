# Plot Method for cv.fdata

For regression: observed vs predicted scatter plot and per-fold residual
boxplots (side by side via patchwork if available, otherwise just the
scatter). For classification: confusion matrix heatmap.

## Usage

``` r
# S3 method for class 'cv.fdata'
plot(x, ...)
```

## Arguments

- x:

  A `cv.fdata` object.

- ...:

  Ignored.

## Value

A `ggplot` object (or a patchwork composition).
