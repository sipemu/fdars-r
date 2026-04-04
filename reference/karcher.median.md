# Karcher Median in Elastic Metric

Compute the Karcher median of functional data in the elastic metric
using the Weiszfeld algorithm. The Karcher median minimises the sum of
distances rather than the sum of squared distances, making it more
robust to outliers than the Karcher mean.

## Usage

``` r
karcher.median(fdataobj, max.iter = 30, tol = 1e-04, lambda = 0, trim = 0)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- max.iter:

  Maximum number of Weiszfeld iterations (default 30).

- tol:

  Convergence tolerance (default 1e-4).

- lambda:

  Regularisation parameter controlling warping smoothness (default 0).

- trim:

  Fraction of most distant curves to trim before aggregation (default 0,
  no trimming).

## Value

An object of class 'karcher.median' with components:

- median:

  fdata of the Karcher median curve

- median_srsf:

  numeric vector of the median SRSF

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

Fletcher, P.T., Venkatasubramanian, S., and Joshi, S. (2009). The
geometric median on Riemannian manifolds with application to robust
atlas estimation. *NeuroImage*, 45(1):S143–S152.

Srivastava, A. and Klassen, E. (2016). *Functional and Shape Data
Analysis*. Springer.

## See also

[`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md),
[`robust.karcher.mean`](https://sipemu.github.io/fdars-r/reference/robust.karcher.mean.md)

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 50)
X <- matrix(0, 15, 50)
for (i in 1:15) X[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
fd <- fdata(X, argvals = t)
km <- karcher.median(fd)
# }
```
