# LIME Explanation

Computes a Local Interpretable Model-agnostic Explanation for a single
observation.

## Usage

``` r
fregre.lime(
  model,
  data,
  observation,
  n.samples = 500,
  kernel.width = 1,
  seed = NULL
)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

- data:

  An fdata object.

- observation:

  Index of the observation to explain (1-based).

- n.samples:

  Number of neighborhood samples (default 500).

- kernel.width:

  Kernel width for local weighting (default 1).

- seed:

  Random seed.

## Value

A list with observation, attributions, local_intercept, local_r_squared,
and kernel_width.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
result <- fregre.lime(fit, fd, observation = 1)
# }
```
