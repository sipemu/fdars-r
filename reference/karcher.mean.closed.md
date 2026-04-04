# Karcher Mean for Closed Curves

Compute the Karcher (Frechet) mean for a sample of closed curves in the
elastic metric. This optimises simultaneously over reparameterisations
and rotations of the parameter domain.

## Usage

``` r
karcher.mean.closed(fdataobj, max.iter = 15, tol = 1e-04, lambda = 0)
```

## Arguments

- fdataobj:

  An object of class 'fdata' containing closed curves.

- max.iter:

  Maximum number of iterations (default 15).

- tol:

  Convergence tolerance (default 1e-4).

- lambda:

  Regularisation parameter controlling warping smoothness (default 0).

## Value

An object of class 'karcher.mean' with components:

- mean:

  fdata of the Karcher mean curve

- mean_srsf:

  numeric vector of the mean SRSF

- aligned:

  fdata of aligned curves

- gammas:

  fdata of warping functions

- n.iter:

  number of iterations used

- converged:

  logical indicating convergence

- fdataobj:

  the original fdata input

## References

Srivastava, A. and Klassen, E. (2016). *Functional and Shape Data
Analysis*. Springer.

## See also

[`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md),
[`elastic.align.pair.closed`](https://sipemu.github.io/fdars-r/reference/elastic.align.pair.closed.md)

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 100)
X <- matrix(0, 10, 100)
for (i in 1:10) X[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
fd <- fdata(X, argvals = t)
km <- karcher.mean.closed(fd)
# }
```
