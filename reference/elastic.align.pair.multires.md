# Multi-Resolution Pairwise Curve Alignment

Align two curves using a coarse-to-fine multi-resolution strategy. The
curves are first coarsened, aligned at low resolution, and then the
alignment is progressively refined to the original grid resolution. This
can improve robustness and speed for highly misaligned curves.

## Usage

``` r
elastic.align.pair.multires(
  f1,
  f2,
  argvals,
  coarsen.factor = 4,
  n.refine = 10,
  step.size = 0.01,
  lambda = 0
)
```

## Arguments

- f1:

  Numeric vector of the first curve's values.

- f2:

  Numeric vector of the second curve's values.

- argvals:

  Numeric vector of the evaluation grid (common to both curves).

- coarsen.factor:

  Coarsening factor for the initial resolution (default 4).

- n.refine:

  Number of refinement steps (default 10).

- step.size:

  Step size for gradient-based refinement (default 0.01).

- lambda:

  Regularisation parameter controlling warping smoothness (default 0).

## Value

A list with components:

- gamma:

  numeric vector of the estimated warping function

- f2.aligned:

  numeric vector of the aligned second curve

- distance:

  elastic distance after alignment

## See also

[`elastic.align`](https://sipemu.github.io/fdars-r/reference/elastic.align.md)

## Examples

``` r
# \donttest{
t <- seq(0, 1, length.out = 200)
f1 <- sin(2 * pi * t)
f2 <- sin(2 * pi * (t^1.5))
res <- elastic.align.pair.multires(f1, f2, t)
# }
```
