# Karcher Mean in Elastic Metric

Compute the Karcher (Fréchet) mean of functional data in the elastic
metric. This simultaneously estimates the mean shape and aligns all
curves.

## Usage

``` r
karcher.mean(
  fdataobj,
  max.iter = 20,
  tol = 1e-04,
  periodic = FALSE,
  rotate.method = "peak",
  rotate.args = list()
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- max.iter:

  Maximum number of iterations (default 20).

- tol:

  Convergence tolerance (default 1e-4).

- periodic:

  Logical; if TRUE, circularly rotate each curve to a canonical position
  before computing the Karcher mean. See
  [`elastic.align`](https://sipemu.github.io/fdars-r/reference/elastic.align.md)
  for details. Default is FALSE.

- rotate.method:

  Rotation method when `periodic = TRUE`: one of `"peak"` (default),
  `"xcorr"`, `"landmark"`, or `"iterative"`. See
  [`periodic.rotate`](https://sipemu.github.io/fdars-r/reference/periodic.rotate.md)
  for details.

- rotate.args:

  A named list of additional arguments passed to
  [`periodic.rotate`](https://sipemu.github.io/fdars-r/reference/periodic.rotate.md)
  (e.g., `reference`, `landmark.func`, `max.iter`).

## Value

An object of class 'karcher.mean' with components:

- mean:

  fdata of the Karcher mean curve

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

- rotations:

  integer vector of circular rotation shifts applied (NULL when periodic
  = FALSE)

- rotate_method:

  the rotation method used (NULL when periodic = FALSE)

## References

Srivastava, A., Klassen, E., Joshi, S.H., and Jermyn, I.H. (2011). Shape
analysis of elastic curves in Euclidean spaces. *IEEE Transactions on
Pattern Analysis and Machine Intelligence*, 33(7):1415–1428.

Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
functional data using phase and amplitude separation. *Computational
Statistics & Data Analysis*, 61:50–66.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
km <- karcher.mean(fd)
# }
```
