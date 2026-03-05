# Alignment Quality Diagnostics

Compute comprehensive alignment quality metrics from a Karcher mean
result, including warp complexity, smoothness, and variance
decomposition.

## Usage

``` r
alignment.quality(fdataobj, karcher)
```

## Arguments

- fdataobj:

  An object of class 'fdata' (the original data).

- karcher:

  An object of class 'karcher.mean'.

## Value

An object of class 'alignment.quality' with components:

- warp_complexity:

  Per-curve geodesic distance from warp to identity

- mean_warp_complexity:

  Mean warp complexity

- warp_smoothness:

  Per-curve bending energy

- mean_warp_smoothness:

  Mean bending energy

- total_variance:

  Total variance of original data

- amplitude_variance:

  Amplitude variance (after alignment)

- phase_variance:

  Phase variance (total - amplitude)

- phase_amplitude_ratio:

  Phase-to-total variance ratio

- pointwise_variance_ratio:

  Pointwise aligned/original variance ratio

- mean_variance_reduction:

  Mean variance reduction

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
km <- karcher.mean(fd, max.iter = 5)
aq <- alignment.quality(fd, km)
# }
```
