# Explanation Stability

Assesses how stable the model explanation is under bootstrap resampling.

## Usage

``` r
fregre.stability(data, y, ncomp, n.boot = 100, seed = NULL)
```

## Arguments

- data:

  An fdata object.

- y:

  Response vector.

- ncomp:

  Number of components.

- n.boot:

  Number of bootstrap resamples (default 100).

- seed:

  Random seed.

## Value

A list with beta_t_std, coefficient_std, metric_std, beta_t_cv,
importance_stability, and n_boot_success.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
result <- fregre.stability(fd, y, ncomp = 3)
# }
```
