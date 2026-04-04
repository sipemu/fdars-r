# Trimmed Karcher Mean in Elastic Metric

Compute a robust version of the Karcher mean by iteratively trimming the
most distant curves before updating the mean estimate. This reduces the
influence of outlying observations on the mean shape.

## Usage

``` r
robust.karcher.mean(
  fdataobj,
  trim = 0.1,
  max.iter = 20,
  tol = 1e-04,
  lambda = 0
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- trim:

  Fraction of most distant curves to trim at each iteration (default
  0.1).

- max.iter:

  Maximum number of iterations (default 20).

- tol:

  Convergence tolerance (default 1e-4).

- lambda:

  Regularisation parameter controlling warping smoothness (default 0).

## Value

An object of class 'karcher.mean' with components:

- mean:

  fdata of the trimmed Karcher mean curve

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

Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
functional data using phase and amplitude separation. *Computational
Statistics & Data Analysis*, 61:50–66.

## See also

[`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md),
[`karcher.median`](https://sipemu.github.io/fdars-r/reference/karcher.median.md)

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 50)
X <- matrix(0, 20, 50)
for (i in 1:20) X[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
fd <- fdata(X, argvals = t)
rkm <- robust.karcher.mean(fd, trim = 0.1)
# }
```
