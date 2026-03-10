# Plot method for peak_detection objects

Plots the original functional data with detected peaks overlaid as red
points.

## Usage

``` r
# S3 method for class 'peak_detection'
plot(x, curves = 1, ...)
```

## Arguments

- x:

  A peak_detection object (from
  [`detect.peaks`](https://sipemu.github.io/fdars-r/reference/detect.peaks.md)).

- curves:

  Index of curve to plot (default: 1).

- ...:

  Additional arguments (ignored).

## Value

A ggplot object.
