# Shape Mean (Karcher Mean in Quotient Space)

Compute the Karcher (Frechet) mean of functional data in a shape
quotient space. This simultaneously estimates the mean shape and aligns
all curves, factoring out the specified nuisance transformations.

## Usage

``` r
shape.mean(
  fdataobj,
  quotient = c("reparameterization", "translation", "scale"),
  lambda = 0,
  max.iter = 20,
  tol = 1e-04
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

- max.iter:

  Maximum number of iterations (default 20).

- tol:

  Convergence tolerance (default 1e-4).

## Value

An object of class `shape.mean` with components:

- mean:

  The shape mean curve (numeric vector)

- mean.srsf:

  SRSF of the mean curve

- gammas:

  Matrix of warping functions (n x m)

- aligned.data:

  Matrix of aligned curves (n x m)

- n.iter:

  Number of iterations

- converged:

  Logical: did the algorithm converge?

- fdataobj:

  Original fdata object

- quotient:

  Quotient space used

- call:

  The matched call

## References

Srivastava, A. and Klassen, E. (2016). *Functional and Shape Data
Analysis*. Springer.

## See also

[`shape.distance`](https://sipemu.github.io/fdars-r/reference/shape.distance.md),
[`shape.representative`](https://sipemu.github.io/fdars-r/reference/shape.representative.md),
[`shape.distance.matrix`](https://sipemu.github.io/fdars-r/reference/shape.distance.matrix.md)

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 50)
X <- matrix(0, 10, 50)
for (i in 1:10) X[i, ] <- sin(2 * pi * (t - i / 50)) + rnorm(50, 0, 0.1)
fd <- fdata(X, argvals = t)
sm <- shape.mean(fd, quotient = "reparameterization", max.iter = 10)
sm
#> Shape Mean (Quotient Space)
#>   Curves: 10 x 50 grid points
#>   Quotient: reparameterization 
#>   Iterations: 10 
#>   Converged: FALSE 
# }
```
