# TSRVF: Transported Square-Root Velocity Function

Functions for mapping aligned functional data to a tangent space on the
Hilbert sphere, enabling standard linear operations (PCA, regression) on
elastically aligned curves. TSRVF Transform

## Usage

``` r
tsrvf.transform(fdataobj, max.iter = 20, tol = 1e-04, lambda = 0)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- max.iter:

  Maximum number of Karcher mean iterations. Default 20.

- tol:

  Convergence tolerance for Karcher mean. Default 1e-4.

- lambda:

  Penalty weight on warp deviation from identity. Default 0.

## Value

An object of class 'tsrvf' with components:

- tangent_vectors:

  fdata of tangent vectors in Euclidean space

- mean:

  fdata of the Karcher mean curve

- gammas:

  fdata of warping functions

- converged:

  logical indicating Karcher mean convergence

- fdataobj:

  the original input

## Details

Compute the TSRVF (Transported SRSF) representation of functional data.
This performs Karcher mean alignment and then maps the aligned SRSFs to
the tangent space at the mean via the inverse exponential map.

## References

Su, J., Kurtek, S., Klassen, E., and Srivastava, A. (2014). Statistical
analysis of trajectories on Riemannian manifolds: bird migration,
hurricane tracking and video surveillance. *Annals of Applied
Statistics*, 8(1):530–552.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
tv <- tsrvf.transform(fd, max.iter = 5)
# }
```
