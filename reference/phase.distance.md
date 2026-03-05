# Phase Distance Matrix

Compute the pairwise phase distance matrix. The phase distance measures
the geodesic distance of the optimal warping function from the identity.

## Usage

``` r
phase.distance(fdataobj, fdataref = NULL, lambda = 0, ...)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- fdataref:

  Optional reference 'fdata'. If NULL, computes self-distances.

- lambda:

  Penalty weight. Default 0.

- ...:

  Additional arguments (ignored).

## Value

A distance matrix.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(60), 6, 10), argvals = seq(0, 1, length.out = 10))
D <- phase.distance(fd)
# }
```
