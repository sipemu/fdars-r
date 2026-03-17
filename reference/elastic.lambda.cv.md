# Cross-Validation for Elastic Alignment Regularization Parameter

Select the optimal regularization parameter (lambda) for elastic
alignment using K-fold cross-validation. The CV error measures how well
the aligned curves generalize to held-out folds.

## Usage

``` r
elastic.lambda.cv(
  fdataobj,
  lambdas = 10^seq(-4, 2, length.out = 20),
  n.folds = 5,
  max.iter = 15,
  tol = 0.001,
  seed = 42
)
```

## Arguments

- fdataobj:

  An object of class `fdata`.

- lambdas:

  Numeric vector of candidate lambda values. Default is
  `10^seq(-4, 2, length.out = 20)`.

- n.folds:

  Number of cross-validation folds (default 5).

- max.iter:

  Maximum iterations for Karcher mean per fold (default 15).

- tol:

  Convergence tolerance (default 1e-3).

- seed:

  Random seed for fold assignment (default 42).

## Value

An object of class `lambda.cv` with components:

- best.lambda:

  The lambda value with the lowest CV error

- cv.scores:

  Numeric vector of CV scores for each lambda

- lambdas:

  The candidate lambda values tested

- call:

  The matched call

## References

Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
functional data using phase and amplitude separation. *Computational
Statistics & Data Analysis*, 61:50–66.

## See also

[`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md),
[`elastic.align`](https://sipemu.github.io/fdars-r/reference/elastic.align.md)

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 50)
X <- matrix(0, 15, 50)
for (i in 1:15) X[i, ] <- sin(2 * pi * (t - i / 60))
fd <- fdata(X, argvals = t)
cv <- elastic.lambda.cv(fd, lambdas = 10^seq(-2, 1, length.out = 10))
cv
#> Lambda Cross-Validation
#>   Candidates tested: 10 
#>   Best lambda: 0.01 
#>   Best CV score: 0.3533 
# }
```
