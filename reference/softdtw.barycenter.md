# Soft-DTW Barycenter

Compute the barycenter (weighted average) of functional data using the
Soft-DTW distance. The barycenter minimizes the sum of Soft-DTW
distances to all input curves.

## Usage

``` r
softdtw.barycenter(fdataobj, gamma = 1, max.iter = 100, tol = 1e-05)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- gamma:

  Smoothing parameter (\> 0). Default is 1.0.

- max.iter:

  Maximum number of iterations. Default is 100.

- tol:

  Convergence tolerance. Default is 1e-5.

## Value

An object of class 'fdata' containing the barycenter curve, with
attributes `n.iter` and `converged`.

## References

Cuturi, M. and Blondel, M. (2017). Soft-DTW: a Differentiable Loss
Function for Time-Series. *Proceedings of the 34th International
Conference on Machine Learning (ICML)*.

## Examples

``` r
fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
bc <- softdtw.barycenter(fd, gamma = 1.0)
```
