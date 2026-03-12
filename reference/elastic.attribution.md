# Elastic PCR Attribution

Decomposes prediction importance into amplitude and phase contributions
after elastic PCR.

## Usage

``` r
elastic.attribution(
  fdataobj,
  y,
  ncomp = 3,
  pca.method = c("vertical", "horizontal", "joint"),
  lambda = 0,
  max.iter = 20,
  tol = 1e-04,
  n.perm = 100,
  seed = NULL
)
```

## Arguments

- fdataobj:

  An fdata object.

- y:

  Response vector.

- ncomp:

  Number of elastic FPCA components.

- pca.method:

  PCA method: "vertical", "horizontal", or "joint".

- lambda:

  Regularization for alignment (default 0).

- max.iter:

  Maximum alignment iterations (default 20).

- tol:

  Convergence tolerance (default 1e-4).

- n.perm:

  Number of permutations (default 100).

- seed:

  Random seed.

## Value

A list with amplitude_r_squared, phase_r_squared, total_r_squared,
amplitude_importance, phase_importance, and p_values.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
result <- elastic.attribution(fd, y, ncomp = 3)
# }
```
