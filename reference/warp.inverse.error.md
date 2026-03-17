# Warp Inverse Error

Compute the maximum reconstruction error of a warp inversion: \\\max_t
\|\gamma(\gamma^{-1}(t)) - t\|\\.

## Usage

``` r
warp.inverse.error(gamma, argvals)
```

## Arguments

- gamma:

  Numeric vector of warping function values.

- argvals:

  Numeric vector of evaluation points.

## Value

A scalar error value.

## See also

[`warp.inverse`](https://sipemu.github.io/fdars-r/reference/warp.inverse.md)

## Examples

``` r
t <- seq(0, 1, length.out = 50)
gamma <- t^2
warp.inverse.error(gamma, t)
#> [1] 1.110223e-16
```
