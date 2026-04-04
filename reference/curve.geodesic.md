# Geodesic Path Between Two Curves

Compute the geodesic (shortest path) between two curves in the elastic
metric. The geodesic is represented as a sequence of intermediate curves
interpolating between the two endpoints in the SRVF space.

## Usage

``` r
curve.geodesic(f1, f2, argvals, n.points = 10, lambda = 0)
```

## Arguments

- f1:

  Numeric vector of the first curve's values.

- f2:

  Numeric vector of the second curve's values.

- argvals:

  Numeric vector of the evaluation grid (common to both curves).

- n.points:

  Number of points along the geodesic to compute (default 10).

- lambda:

  Regularisation parameter controlling warping smoothness (default 0).

## Value

An object of class 'curve.geodesic' with components:

- path:

  fdata of the interpolated curves along the geodesic

- distance:

  the total geodesic (elastic) distance between f1 and f2

- gamma:

  numeric vector of the optimal warping function

- argvals:

  the evaluation grid

## References

Srivastava, A. and Klassen, E. (2016). *Functional and Shape Data
Analysis*. Springer.

## See also

[`elastic.distance`](https://sipemu.github.io/fdars-r/reference/elastic.distance.md)

## Examples

``` r
# \donttest{
t <- seq(0, 1, length.out = 100)
f1 <- sin(2 * pi * t)
f2 <- cos(2 * pi * t)
geo <- curve.geodesic(f1, f2, t, n.points = 8)
# }
```
