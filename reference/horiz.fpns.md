# Horizontal Functional Principal Nested Spheres

Perform horizontal functional PCA using the principal nested spheres
(fPNS) framework on warping functions. This captures the main modes of
phase variability from the warping functions estimated during Karcher
mean computation.

## Usage

``` r
horiz.fpns(karcher, n.components = 3)
```

## Arguments

- karcher:

  An object of class 'karcher.mean' (result of
  [`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md)).

- n.components:

  Number of principal components to retain (default 3).

## Value

A list with components:

- scores:

  matrix of phase component scores (curves x components)

- eigenvalues:

  numeric vector of eigenvalues

- eigenfunctions:

  matrix of phase eigenfunctions (components x grid)

- variance.explained:

  numeric vector of proportion of variance explained by each component

## References

Jung, S., Dryden, I.L., and Marron, J.S. (2012). Analysis of principal
nested spheres. *Biometrika*, 99(3):551–568.

Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
functional data using phase and amplitude separation. *Computational
Statistics & Data Analysis*, 61:50–66.

## See also

[`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md),
[`gauss.model`](https://sipemu.github.io/fdars-r/reference/gauss.model.md)

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 50)
X <- matrix(0, 20, 50)
for (i in 1:20) X[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
fd <- fdata(X, argvals = t)
km <- karcher.mean(fd)
hfp <- horiz.fpns(km, n.components = 2)
# }
```
