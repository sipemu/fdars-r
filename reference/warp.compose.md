# Compose Two Warping Functions

Computes the composition (gamma1 o gamma2)(t) = gamma1(gamma2(t)).

## Usage

``` r
warp.compose(gamma1, gamma2, fdataobj)
```

## Arguments

- gamma1:

  Numeric vector — outer warping function.

- gamma2:

  Numeric vector — inner warping function.

- fdataobj:

  An fdata object providing the grid (argvals).

## Value

Numeric vector of the composed warping function.

## Examples

``` r
# \donttest{
t <- seq(0, 1, length.out = 20)
gamma1 <- t^2
gamma2 <- sqrt(t)
composed <- warp.compose(gamma1, gamma2, t)
# }
```
