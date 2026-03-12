# FPC SHAP Values

Computes Shapley values for the FPC scores.

## Usage

``` r
fregre.shap(model, data, n.samples = 100, seed = NULL)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- data:

  An fdata object (the training data).

- n.samples:

  Number of samples for logistic SHAP (default 100).

- seed:

  Random seed (used for logistic models).

## Value

A list with values (SHAP matrix), base_value, and mean_scores.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
result <- fregre.shap(fit, fd)
# }
```
