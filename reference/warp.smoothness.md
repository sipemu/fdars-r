# Warping Smoothness

Computes the bending energy of a warping function, measuring its
roughness.

## Usage

``` r
warp.smoothness(gamma, argvals)
```

## Arguments

- gamma:

  Numeric vector of warping function values.

- argvals:

  Numeric vector of grid points.

## Value

A scalar smoothness (bending energy) value.

## Examples

``` r
t <- seq(0, 1, length.out = 20)
warp.smoothness(t^2, t)
#> [1] 4
```
