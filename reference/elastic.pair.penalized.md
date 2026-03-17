# Penalized Elastic Pairwise Alignment

Aligns curve `f2` to curve `f1` using elastic alignment with a
configurable roughness penalty on the warping function.

## Usage

``` r
elastic.pair.penalized(
  f1,
  f2,
  argvals = NULL,
  penalty = c("first_order", "second_order", "combined"),
  lambda = 1
)
```

## Arguments

- f1:

  Numeric vector (reference curve).

- f2:

  Numeric vector (curve to align).

- argvals:

  Numeric vector of evaluation points. If `NULL`, defaults to
  `seq(0, 1, length.out = length(f1))`.

- penalty:

  Character: penalty type. One of `"first_order"`, `"second_order"`, or
  `"combined"`.

- lambda:

  Numeric penalty weight (default 1.0).

## Value

A list with components:

- gamma:

  Warping function

- f.aligned:

  Aligned version of f2

- distance:

  Elastic distance after alignment

## See also

[`elastic.pair`](https://sipemu.github.io/fdars-r/reference/elastic.pair.md),
[`elastic.lambda.cv`](https://sipemu.github.io/fdars-r/reference/elastic.lambda.cv.md)

## Examples

``` r
# \donttest{
t <- seq(0, 1, length.out = 50)
f1 <- sin(2 * pi * t)
f2 <- sin(2 * pi * t^1.5)
res <- elastic.pair.penalized(f1, f2, argvals = t)
# }
```
