# Elastic Logistic Classification

Fits a functional logistic regression with elastic alignment.

## Usage

``` r
elastic.logistic(
  fdataobj,
  y,
  ncomp.beta = 10,
  lambda = 0,
  max.iter = 20,
  tol = 1e-04
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- y:

  Binary response vector (0/1).

- ncomp.beta:

  Number of basis functions for beta (default 10).

- lambda:

  Regularization parameter (default 0).

- max.iter:

  Maximum iterations (default 20).

- tol:

  Convergence tolerance (default 1e-4).

## Value

An object of class 'elastic.logistic' with components:

- alpha:

  Intercept

- beta:

  Beta coefficient function

- probabilities:

  Predicted probabilities

- predicted.classes:

  Predicted class labels (0/1)

- accuracy:

  Classification accuracy

- loss:

  Final loss value

- gammas:

  Estimated warping functions

- aligned.srsfs:

  Aligned SRSF transforms

- n.iter:

  Number of iterations

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- sample(0:1, 50, replace = TRUE)
fit <- elastic.logistic(fd, y)
fit
#> Elastic Logistic Classification
#>   Accuracy: 88 %
#>   Iterations: 20 
# }
```
