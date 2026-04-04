# Elastic Distance Between Closed Curves

Compute the elastic (Fisher-Rao) distance between two closed curves,
optimising over both reparameterisation and rotation.

## Usage

``` r
elastic.distance.closed(f1, f2, argvals, lambda = 0)
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

Numeric scalar; the elastic distance between the two closed curves.

## References

Srivastava, A. and Klassen, E. (2016). *Functional and Shape Data
Analysis*. Springer.

## See also

[`elastic.align.pair.closed`](https://sipemu.github.io/fdars-r/reference/elastic.align.pair.closed.md),
[`karcher.mean.closed`](https://sipemu.github.io/fdars-r/reference/karcher.mean.closed.md)

## Examples

``` r
# \donttest{
t <- seq(0, 1, length.out = 100)
f1 <- sin(2 * pi * t)
f2 <- sin(2 * pi * t + 0.5)
d <- elastic.distance.closed(f1, f2, t)
# }
```
