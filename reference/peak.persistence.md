# Peak Persistence Analysis

Analyse the persistence of peaks (local maxima) across a range of
regularisation strengths. As the regularisation parameter \\\lambda\\
increases, aligned curves become smoother and spurious peaks disappear.
This function tracks which peaks persist across the regularisation path,
providing a multi-scale view of the data's peak structure.

## Usage

``` r
peak.persistence(fdataobj, lambdas = NULL, max.iter = 10, tol = 1e-04)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- lambdas:

  Numeric vector of regularisation values to sweep. If NULL (default),
  uses `seq(0, 5, length.out = 30)`.

- max.iter:

  Maximum number of Karcher mean iterations at each lambda (default 10).

- tol:

  Convergence tolerance (default 1e-4).

## Value

A list with components:

- lambdas:

  the regularisation values used

- n_peaks:

  integer vector of peak counts at each lambda

- peak_locations:

  list of numeric vectors giving peak locations at each lambda

- persistence:

  numeric vector of persistence scores for each peak

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
t <- seq(0, 1, length.out = 100)
X <- matrix(0, 15, 100)
for (i in 1:15) X[i, ] <- sin(4 * pi * t + runif(1, -0.3, 0.3))
fd <- fdata(X, argvals = t)
pp <- peak.persistence(fd)
# }
```
