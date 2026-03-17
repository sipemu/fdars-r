# Invert a Warping Function

Compute the functional inverse of a warping function \\\gamma\\, i.e.,
find \\\gamma^{-1}\\ such that \\\gamma(\gamma^{-1}(t)) = t\\.

## Usage

``` r
warp.inverse(gamma, argvals)
```

## Arguments

- gamma:

  Numeric vector of warping function values.

- argvals:

  Numeric vector of evaluation points.

## Value

Numeric vector of the inverse warping function.

## See also

[`warp.inverse.error`](https://sipemu.github.io/fdars-r/reference/warp.inverse.error.md),
[`warp.compose`](https://sipemu.github.io/fdars-r/reference/warp.compose.md)

## Examples

``` r
t <- seq(0, 1, length.out = 50)
gamma <- t^2
gamma_inv <- warp.inverse(gamma, t)
```
