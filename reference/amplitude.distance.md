# Amplitude Distance Matrix

Compute the pairwise amplitude distance matrix. The amplitude distance
is the elastic (Fisher-Rao) distance after optimal alignment.

## Usage

``` r
amplitude.distance(fdataobj, fdataref = NULL, lambda = 0, ...)
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
D <- amplitude.distance(fd)
# }
```
