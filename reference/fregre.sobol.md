# Sobol Indices

Computes first-order and total-order Sobol sensitivity indices for each
FPC component.

## Usage

``` r
fregre.sobol(model, data, y = NULL, n.samples = 500, seed = NULL)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- data:

  An fdata object (the training data).

- y:

  Response vector (used for lm models).

- n.samples:

  Number of samples for logistic Sobol (default 500).

- seed:

  Random seed (used for logistic models).

## Value

A list with first_order, total_order, var_y, and component_variance.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
result <- fregre.sobol(fit, fd, y)
# }
```
