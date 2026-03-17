# Pairwise Alignment Diagnostics

Compute diagnostic information for a single pairwise alignment result.

## Usage

``` r
alignment.diagnostics.pairwise(f1, f2, result, argvals)
```

## Arguments

- f1:

  Numeric vector (reference curve).

- f2:

  Numeric vector (curve that was aligned).

- result:

  A pairwise alignment result list with components `$gamma` (or
  `$gamma`), `$f.aligned` (or `$f_aligned`), and `$distance`.

- argvals:

  Numeric vector of evaluation points.

## Value

A list with diagnostic fields: `curve.index`, `warp.complexity`,
`warp.smoothness`, `is.under.aligned`, `is.over.aligned`,
`has.non.monotone`, `residual`, `distance.ratio`, `flagged`, `issues`.

## See also

[`elastic.pair`](https://sipemu.github.io/fdars-r/reference/elastic.pair.md),
[`alignment.diagnostics`](https://sipemu.github.io/fdars-r/reference/alignment.diagnostics.md)

## Examples

``` r
# \donttest{
t <- seq(0, 1, length.out = 50)
f1 <- sin(2 * pi * t)
f2 <- sin(2 * pi * t^1.5)
fd1 <- fdata(matrix(f1, 1, 50), argvals = t)
fd2 <- fdata(matrix(f2, 1, 50), argvals = t)
res <- elastic.pair(fd1, fd2)
diag <- alignment.diagnostics.pairwise(f1, f2, res, t)
# }
```
