# LRT Outlier Threshold with Bootstrap Distribution

Computes the bootstrap threshold AND the full sorted null distribution
for LRT-based outlier detection. The returned distribution enables
per-curve p-value computation:
`p = (sum(boot_dist >= d) + 1) / (B + 1)`.

## Usage

``` r
outliers.lrt.dist(
  fdataobj,
  nb = 200,
  smo = 0.05,
  trim = 0.1,
  seed = NULL,
  percentile = 0.99
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- nb:

  Number of bootstrap replications (default 200).

- smo:

  Smoothing parameter for bootstrap noise (default 0.05).

- trim:

  Proportion of curves to trim (default 0.1).

- seed:

  Random seed for reproducibility.

- percentile:

  Percentile for threshold (default 0.99).

## Value

A list with components:

- threshold:

  Threshold at the specified percentile

- boot.distribution:

  Sorted bootstrap null distribution of max-distances

## Examples

``` r
t <- seq(0, 1, length.out = 50)
X <- matrix(0, 30, 50)
for (i in 1:30) X[i, ] <- sin(2*pi*t) + rnorm(50, sd = 0.1)
fd <- fdata(X, argvals = t)
res <- outliers.lrt.dist(fd, nb = 100)
res$threshold
#> [1] 1.929467
```
