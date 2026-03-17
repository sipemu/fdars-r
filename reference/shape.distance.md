# Shape Distance Between Two Curves

Compute the elastic shape distance between two curves in a specified
quotient space. The shape distance factors out the specified nuisance
transformations (reparameterization, translation, scale).

## Usage

``` r
shape.distance(
  f1,
  f2,
  argvals = NULL,
  quotient = c("reparameterization", "translation", "scale"),
  lambda = 0
)
```

## Arguments

- f1:

  Numeric vector (first curve).

- f2:

  Numeric vector (second curve).

- argvals:

  Numeric vector of evaluation points. If `NULL`, defaults to
  `seq(0, 1, length.out = length(f1))`.

- quotient:

  Character: the quotient space. One of `"reparameterization"`,
  `"translation"`, or `"scale"`.

- lambda:

  Regularization parameter (default 0).

## Value

A list with components:

- distance:

  Shape distance

- gamma:

  Optimal warping function

- f2.aligned:

  Aligned version of f2

## References

Srivastava, A. and Klassen, E. (2016). *Functional and Shape Data
Analysis*. Springer.

## See also

[`shape.representative`](https://sipemu.github.io/fdars-r/reference/shape.representative.md),
[`shape.mean`](https://sipemu.github.io/fdars-r/reference/shape.mean.md),
[`shape.distance.matrix`](https://sipemu.github.io/fdars-r/reference/shape.distance.matrix.md)

## Examples

``` r
# \donttest{
t <- seq(0, 1, length.out = 50)
f1 <- sin(2 * pi * t)
f2 <- 2 * sin(2 * pi * t^1.5) + 1
d <- shape.distance(f1, f2, argvals = t, quotient = "scale")
d$distance
#> [1] 0.2732183
# }
```
