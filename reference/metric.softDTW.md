# Soft Dynamic Time Warping Distance

Computes the Soft-DTW distance matrix between functional data objects.
Soft-DTW is a differentiable relaxation of DTW that uses a smoothing
parameter gamma to control the softness of the minimum operation.

## Usage

``` r
metric.softDTW(fdataobj, fdataref = NULL, gamma = 1, divergence = FALSE, ...)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- fdataref:

  An object of class 'fdata'. If NULL, computes self-distances.

- gamma:

  Smoothing parameter (\> 0). Smaller values approximate hard DTW.
  Default is 1.0.

- divergence:

  Logical. If TRUE, computes the Soft-DTW divergence (non-negative, zero
  for identical series) instead of raw Soft-DTW distance. Default is
  FALSE.

- ...:

  Additional arguments (ignored).

## Value

A distance matrix.

## References

Cuturi, M. and Blondel, M. (2017). Soft-DTW: a Differentiable Loss
Function for Time-Series. *Proceedings of the 34th International
Conference on Machine Learning (ICML)*.

## Examples

``` r
fd <- fdata(matrix(rnorm(100), 10, 10))
D <- metric.softDTW(fd, gamma = 1.0)
D_div <- metric.softDTW(fd, gamma = 1.0, divergence = TRUE)
```
