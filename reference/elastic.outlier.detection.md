# Elastic Outlier Detection

Identify outlying curves using SRVF-based elastic distances and Tukey's
fence rule. Outliers are detected in both the amplitude and phase
domains.

## Usage

``` r
elastic.outlier.detection(
  fdataobj,
  lambda = 0,
  alpha = 0.05,
  use.median = FALSE
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- lambda:

  Regularisation parameter controlling warping smoothness (default 0).

- alpha:

  Significance level for the Tukey fence (default 0.05).

- use.median:

  Logical; if TRUE, use the Karcher median instead of the Karcher mean
  as the central tendency reference (default FALSE).

## Value

A list with components:

- outlier_indices:

  integer vector of detected outlier indices

- amplitude_outliers:

  integer vector of amplitude outlier indices

- phase_outliers:

  integer vector of phase outlier indices

- amplitude_distances:

  numeric vector of amplitude distances

- phase_distances:

  numeric vector of phase distances

- amplitude_fence:

  numeric; the amplitude fence threshold

- phase_fence:

  numeric; the phase fence threshold

## References

Xie, W., Kurtek, S., Bharath, K., and Sun, Y. (2017). A geometric
approach to visualization of variability in functional data. *Journal of
the American Statistical Association*, 112(519):979–993.

## See also

[`elastic.depth`](https://sipemu.github.io/fdars-r/reference/elastic.depth.md)

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 50)
X <- matrix(0, 20, 50)
for (i in 1:20) X[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
# Add an outlier
X[20, ] <- 5 * cos(4 * pi * t)
fd <- fdata(X, argvals = t)
out <- elastic.outlier.detection(fd)
# }
```
