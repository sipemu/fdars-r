# Shape Distance Matrix

Compute the pairwise shape distance matrix for functional data in a
specified quotient space.

## Usage

``` r
shape.distance.matrix(
  fdataobj,
  quotient = c("reparameterization", "translation", "scale"),
  lambda = 0
)
```

## Arguments

- fdataobj:

  An object of class `fdata`.

- quotient:

  Character: the quotient space. One of `"reparameterization"`,
  `"translation"`, or `"scale"`.

- lambda:

  Regularization parameter (default 0).

## Value

A symmetric distance matrix (n x n).

## References

Srivastava, A. and Klassen, E. (2016). *Functional and Shape Data
Analysis*. Springer.

## See also

[`shape.distance`](https://sipemu.github.io/fdars-r/reference/shape.distance.md),
[`shape.mean`](https://sipemu.github.io/fdars-r/reference/shape.mean.md),
[`elastic.distance`](https://sipemu.github.io/fdars-r/reference/elastic.distance.md)

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 30)
X <- matrix(0, 8, 30)
for (i in 1:8) X[i, ] <- sin(2 * pi * (t - i / 40))
fd <- fdata(X, argvals = t)
D <- shape.distance.matrix(fd, quotient = "reparameterization")
# }
```
