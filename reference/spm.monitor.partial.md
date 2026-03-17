# Monitor a Partially-Observed Curve

Monitors a single curve that is only partially observed on the domain.
The unobserved portion is completed using one of several methods before
projecting onto the Phase I FPCA basis.

## Usage

``` r
spm.monitor.partial(
  chart,
  partial.data,
  n.observed,
  completion = c("conditional", "projection", "zero_pad"),
  ncomp = NULL,
  alpha = 0.05
)
```

## Arguments

- chart:

  An object of class `spm.chart` from
  [`spm.phase1`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md).

- partial.data:

  Numeric vector of observed values for one curve.

- n.observed:

  Integer; how many of the grid points have been observed (from the
  start of the domain).

- completion:

  Completion method: `"conditional"` (conditional expectation / BLUP),
  `"projection"` (partial projection), or `"zero_pad"`.

- ncomp:

  Number of components (default: same as chart).

- alpha:

  Significance level (default 0.05).

## Value

A list with components:

- scores:

  Estimated FPC scores

- t2:

  T-squared statistic

- t2.alarm:

  Logical: TRUE if T-squared exceeds UCL

- domain.fraction:

  Fraction of domain observed

- completed.curve:

  The completed curve (full domain), or NULL

## See also

[`spm.monitor`](https://sipemu.github.io/fdars-r/reference/spm.monitor.md)
for fully-observed curves,
[`spm.monitor.partial.batch`](https://sipemu.github.io/fdars-r/reference/spm.monitor.partial.batch.md)
for multiple partial curves

## Examples

``` r
# \donttest{
set.seed(1)
n <- 50; m <- 30
argvals <- seq(0, 1, length.out = m)
X <- matrix(rnorm(n * m), n, m)
fd <- fdata(X, argvals = argvals)
chart <- spm.phase1(fd, ncomp = 3)

# Monitor a partially observed curve (first 20 of 30 grid points)
partial <- rnorm(30)
result <- spm.monitor.partial(chart, partial, n.observed = 20)
result$t2
#> [1] 1.220404
result$domain.fraction
#> [1] 0.6551724
# }
```
