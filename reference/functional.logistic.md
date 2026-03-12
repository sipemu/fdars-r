# Functional Logistic Regression

Binary logistic regression with functional and optional scalar
predictors using FPC projection and IRLS.

## Usage

``` r
functional.logistic(
  fdataobj,
  y,
  scalar.covariates = NULL,
  ncomp = 3,
  max.iter = 100,
  tol = 1e-06
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- y:

  Binary response vector (0/1).

- scalar.covariates:

  Optional matrix of scalar covariates.

- ncomp:

  Number of FPC components (default 3).

- max.iter:

  Maximum IRLS iterations (default 100).

- tol:

  Convergence tolerance (default 1e-6).

## Value

A fitted object of class 'fregre.logistic' with components:

- probabilities:

  Predicted probabilities P(Y=1)

- predicted.classes:

  Predicted class labels (0 or 1)

- accuracy:

  Classification accuracy on training data

- log.likelihood:

  Log-likelihood at convergence

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- factor(sample(0:1, 50, replace = TRUE))
result <- functional.logistic(fd, y, ncomp = 3)
#> Error in functional.logistic(fd, y, ncomp = 3): functional.logistic failed: check data dimensions and response
result$accuracy
#> Error: object 'result' not found
# }
```
