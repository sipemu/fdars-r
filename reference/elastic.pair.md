# Elastic Pairwise Alignment

Aligns curve f2 to curve f1 using elastic (Fisher-Rao) alignment.
Returns the aligned curve, warping function, and elastic distance.

## Usage

``` r
elastic.pair(f1, f2)
```

## Arguments

- f1:

  An fdata object (reference curve, first row used).

- f2:

  An fdata object (curve to align, first row used).

## Value

A list with components:

- `f.aligned` — fdata of the aligned curve.

- `gamma` — Numeric warping function.

- `distance` — Elastic distance after alignment.

## Examples

``` r
# \donttest{
t <- seq(0, 1, length.out = 50)
f1 <- fdata(matrix(sin(2 * pi * t), 1, 50), argvals = t)
f2 <- fdata(matrix(sin(2 * pi * t^1.5), 1, 50), argvals = t)
res <- elastic.pair(f1, f2)
# }
```
