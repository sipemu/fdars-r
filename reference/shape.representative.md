# Shape Representative (Orbit Representative)

Compute the canonical representative of a curve in a quotient space.
Depending on the quotient, this factors out reparameterization,
translation, and/or scale.

## Usage

``` r
shape.representative(
  f,
  argvals = NULL,
  quotient = c("reparameterization", "translation", "scale")
)
```

## Arguments

- f:

  Numeric vector of curve values.

- argvals:

  Numeric vector of evaluation points. If `NULL`, defaults to
  `seq(0, 1, length.out = length(f))`.

- quotient:

  Character: the quotient space. One of `"reparameterization"`,
  `"translation"`, or `"scale"`.

## Value

A list with components:

- representative:

  The canonical representative curve

- representative.srsf:

  SRSF of the representative

- gamma:

  Warping function to the representative

- translation:

  Translation applied (0 for reparameterization-only)

- scale:

  Scale factor applied (1 for reparameterization-only)

## References

Srivastava, A. and Klassen, E. (2016). *Functional and Shape Data
Analysis*. Springer.

## See also

[`shape.distance`](https://sipemu.github.io/fdars-r/reference/shape.distance.md),
[`shape.mean`](https://sipemu.github.io/fdars-r/reference/shape.mean.md)

## Examples

``` r
# \donttest{
t <- seq(0, 1, length.out = 50)
f <- sin(2 * pi * t)
rep <- shape.representative(f, argvals = t)
# }
```
