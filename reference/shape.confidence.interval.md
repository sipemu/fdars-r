# Bootstrap Confidence Bands for Elastic Karcher Mean

Construct pointwise bootstrap confidence bands for the elastic Karcher
mean by resampling curves and recomputing the mean for each bootstrap
replicate.

## Usage

``` r
shape.confidence.interval(
  fdataobj,
  n.bootstrap = 200,
  confidence.level = 0.95,
  lambda = 0,
  max.iter = 15,
  tol = 1e-04,
  seed = 42
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- n.bootstrap:

  Number of bootstrap replicates (default 200).

- confidence.level:

  Confidence level for the bands (default 0.95).

- lambda:

  Regularisation parameter controlling warping smoothness (default 0).

- max.iter:

  Maximum number of Karcher mean iterations per replicate (default 15).

- tol:

  Convergence tolerance for each Karcher mean computation (default
  1e-4).

- seed:

  Random seed for reproducibility (default 42).

## Value

An object of class 'shape.ci' with components:

- mean:

  fdata of the Karcher mean curve

- lower:

  numeric vector of lower confidence band values

- upper:

  numeric vector of upper confidence band values

- confidence.level:

  the confidence level used

- n.bootstrap:

  the number of bootstrap replicates used

- argvals:

  the evaluation grid

## References

Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
functional data using phase and amplitude separation. *Computational
Statistics & Data Analysis*, 61:50–66.

## See also

[`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md)

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 50)
X <- matrix(0, 20, 50)
for (i in 1:20) X[i, ] <- sin(2 * pi * t + runif(1, -0.2, 0.2))
fd <- fdata(X, argvals = t)
ci <- shape.confidence.interval(fd, n.bootstrap = 50)
# }
```
