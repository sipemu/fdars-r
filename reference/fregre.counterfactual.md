# Counterfactual Explanation

Finds the minimal change in FPC scores needed to achieve a target
prediction value. For logistic models, uses gradient-based optimization
to find the class-flipping counterfactual.

## Usage

``` r
fregre.counterfactual(
  model,
  data,
  observation,
  target.value = NULL,
  max.iter = 100,
  step.size = 0.1
)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- data:

  An fdata object.

- observation:

  Index of the observation (1-based).

- target.value:

  Target prediction value (for lm models).

- max.iter:

  Maximum iterations for logistic counterfactual (default 100).

- step.size:

  Step size for logistic counterfactual (default 0.1).

## Value

A list with observation, original_scores, counterfactual_scores,
delta_scores, delta_function, distance, original_prediction,
counterfactual_prediction, and found.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
result <- fregre.counterfactual(fit, fd, observation = 1, target.value = 0)
# }
```
