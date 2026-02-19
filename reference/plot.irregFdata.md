# Plot method for irregFdata objects

Plot method for irregFdata objects

## Usage

``` r
# S3 method for class 'irregFdata'
plot(
  x,
  ...,
  col = NULL,
  lty = 1,
  lwd = 1,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  add = FALSE,
  alpha = 0.7
)
```

## Arguments

- x:

  An irregFdata object.

- ...:

  Additional arguments passed to `plot`.

- col:

  Colors for curves.

- lty:

  Line type.

- lwd:

  Line width.

- main:

  Plot title.

- xlab:

  X-axis label.

- ylab:

  Y-axis label.

- add:

  Logical. If TRUE, add to existing plot.

- alpha:

  Transparency for many curves.

## Value

Invisibly returns the input object `x`.
