# Landmark-Constrained Elastic Alignment

Align a set of curves to a target curve with landmark constraints. The
alignment is performed using segmented dynamic programming that passes
through specified landmark positions.

## Usage

``` r
elastic.align.constrained(
  fdataobj,
  target = NULL,
  landmark.pairs = NULL,
  kind = "peak",
  min.prominence = 0,
  expected.count = 0,
  lambda = 0
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- target:

  Target curve (numeric vector or single-curve fdata). If NULL, uses the
  cross-sectional mean.

- landmark.pairs:

  A two-column matrix of landmark pairs (target_t, source_t) for
  explicit constraints. If NULL, auto-detects.

- kind:

  Landmark type for auto-detection: "peak", "valley", "zero",
  "inflection". Default "peak".

- min.prominence:

  Minimum prominence for auto-detection. Default 0.

- expected.count:

  Expected number of landmarks for auto-detection. Default 0 (all
  detected).

- lambda:

  Penalty weight on warp deviation from identity. Default 0.

## Value

An object of class 'elastic.align' with components:

- aligned:

  fdata of aligned curves

- gammas:

  fdata of warping functions

- distances:

  numeric vector of elastic distances

- target:

  the target curve

- fdataobj:

  the original input

## Examples

``` r
# \donttest{
t <- seq(0, 1, length.out = 100)
X <- matrix(0, 10, 100)
for (i in 1:10) X[i, ] <- sin(2*pi*(t - i/50))
fd <- fdata(X, argvals = t)
res <- elastic.align.constrained(fd, kind = "peak")
# }
```
