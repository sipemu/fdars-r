# Alignment Diagnostics

Compute per-curve diagnostic information from a Karcher mean alignment,
identifying curves that may be poorly aligned (under-aligned,
over-aligned, or with non-monotone warps).

## Usage

``` r
alignment.diagnostics(fdataobj, karcher.result)
```

## Arguments

- fdataobj:

  An object of class `fdata` (original data).

- karcher.result:

  An object of class `karcher.mean`.

## Value

An object of class `alignment.diagnostics` with components:

- diagnostics:

  List of per-curve diagnostic records, each containing `curve.index`,
  `warp.complexity`, `warp.smoothness`, `is.under.aligned`,
  `is.over.aligned`, `has.non.monotone`, `residual`, `distance.ratio`,
  `flagged`, and `issues`

- flagged.indices:

  Integer vector of 1-indexed flagged curve indices

- n.flagged:

  Number of flagged curves

- health.score:

  Overall alignment health score in \\\[0, 1\]\\

- call:

  The matched call

## See also

[`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md),
[`alignment.quality`](https://sipemu.github.io/fdars-r/reference/alignment.quality.md),
[`alignment.diagnostics.pairwise`](https://sipemu.github.io/fdars-r/reference/alignment.diagnostics.pairwise.md)

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 50)
X <- matrix(0, 15, 50)
for (i in 1:15) X[i, ] <- sin(2 * pi * (t - i / 60))
fd <- fdata(X, argvals = t)
km <- karcher.mean(fd, max.iter = 5)
diag <- alignment.diagnostics(fd, km)
diag
#> Alignment Diagnostics
#>   Curves: 15 
#>   Flagged: 4 of 15 (26.7%) 
#>   Health score: 0.7333 
#>   Flagged indices: 8, 9, 10, 14 
# }
```
