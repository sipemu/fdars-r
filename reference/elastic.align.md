# Elastic Curve Alignment

Align a set of functional curves using the elastic (Fisher-Rao)
framework. When no target is specified, the Karcher mean is used as the
alignment target.

## Usage

``` r
elastic.align(fdataobj, target = NULL, periodic = FALSE)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- target:

  Optional target curve (numeric vector). If NULL, uses the
  cross-sectional mean as the alignment target.

- periodic:

  Logical; if TRUE, circularly rotate each curve to a canonical position
  before elastic alignment. This two-stage approach handles periodic
  functional data (e.g., data on \\\[0, 2\pi\]\\ where \\f(0) =
  f(2\pi)\\) that would otherwise be poorly aligned due to fixed
  boundary constraints \\\gamma(0)=0, \gamma(1)=1\\. Default is FALSE.

## Value

An object of class 'elastic.align' with components:

- aligned:

  fdata of aligned curves

- gammas:

  fdata of warping functions

- distances:

  numeric vector of elastic distances

- target:

  the target curve used for alignment

- fdataobj:

  the original fdata input

- rotations:

  integer vector of circular rotation shifts applied (NULL when periodic
  = FALSE)

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
res <- elastic.align(fd)
# }
```
