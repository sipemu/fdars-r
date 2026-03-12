# Warping Complexity

Computes the geodesic distance of a warping function from the identity,
measuring how much phase variability the warp introduces.

## Usage

``` r
warp.complexity(gamma, argvals)
```

## Arguments

- gamma:

  Numeric vector of warping function values.

- argvals:

  Numeric vector of grid points.

## Value

A scalar complexity value (0 = identity warp).

## Examples

``` r
t <- seq(0, 1, length.out = 20)
warp.complexity(t^2, t)
#> [1] 0.349855
warp.complexity(t, t)  # identity = 0
#> [1] 0
```
