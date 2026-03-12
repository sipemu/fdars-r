# Anchor Explanation

Finds a rule (set of FPC score conditions) that anchors the prediction
to a narrow range with high precision.

## Usage

``` r
fregre.anchor(model, data, observation, precision = 0.95, n.bins = 10)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- data:

  An fdata object.

- observation:

  Index of the observation (1-based).

- precision:

  Minimum precision threshold (default 0.95).

- n.bins:

  Number of bins for discretizing scores (default 10).

## Value

A list with conditions, precision, coverage, n_matching, observation,
and predicted_value.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
result <- fregre.anchor(fit, fd, observation = 1)
# }
```
