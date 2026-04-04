# Elastic Alignment of Closed Curves

Align two closed curves in the elastic framework. Closed curves satisfy
\\f(0) = f(1)\\ and the alignment optimises over both reparameterisation
and rotation (circular shift) of the parameter domain.

## Usage

``` r
elastic.align.pair.closed(f1, f2, argvals, lambda = 0)
```

## Arguments

- f1:

  Numeric vector of the first closed curve's values.

- f2:

  Numeric vector of the second closed curve's values.

- argvals:

  Numeric vector of the evaluation grid (common to both curves).

- lambda:

  Regularisation parameter controlling warping smoothness (default 0).

## Value

A list with components:

- gamma:

  numeric vector of the estimated warping function

- f2.aligned:

  numeric vector of the aligned second curve

- distance:

  elastic distance after alignment

- rotation:

  the optimal rotation applied

## References

Srivastava, A. and Klassen, E. (2016). *Functional and Shape Data
Analysis*. Springer.

## See also

[`elastic.distance.closed`](https://sipemu.github.io/fdars-r/reference/elastic.distance.closed.md),
[`karcher.mean.closed`](https://sipemu.github.io/fdars-r/reference/karcher.mean.closed.md)

## Examples

``` r
# \donttest{
t <- seq(0, 1, length.out = 100)
f1 <- sin(2 * pi * t)
f2 <- sin(2 * pi * t + 0.5)
res <- elastic.align.pair.closed(f1, f2, t)
# }
```
