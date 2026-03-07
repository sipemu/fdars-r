# Periodic Rotation for Functional Data

Circularly rotate each curve to a canonical position before elastic
alignment of periodic functional data (e.g., data on \\\[0, 2\pi\]\\
where \\f(0) = f(2\pi)\\).

## Usage

``` r
periodic.rotate(
  fdataobj,
  target.pos = NULL,
  method = c("peak", "xcorr", "landmark", "iterative"),
  reference = NULL,
  landmark.func = NULL,
  max.iter = 5,
  ...
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- target.pos:

  Target grid index for the feature position. If NULL (default), uses
  `floor(ncol(fdataobj$data) / 4)`.

- method:

  Rotation method: `"peak"` (default), `"xcorr"`, `"landmark"`, or
  `"iterative"`.

- reference:

  Reference curve for cross-correlation methods. A numeric vector or
  single-curve `fdata`. If NULL, uses the cross-sectional mean.

- landmark.func:

  A function taking a numeric vector (one curve) and returning a single
  integer grid index. Required when `method = "landmark"`.

- max.iter:

  Maximum number of rotation-alignment iterations for
  `method = "iterative"`. Default 5.

- ...:

  Additional arguments (ignored).

## Value

A list with components:

- fdataobj:

  fdata of rotated curves

- shifts:

  integer vector of circular shifts applied to each curve

## Details

Four methods are available:

- `"peak"`:

  Shift each curve's global maximum to `target.pos` using `which.max`.
  Fast and simple, but fragile when curves have multiple peaks of
  similar height.

- `"xcorr"`:

  Circular cross-correlation with a reference curve (via
  [`stats::convolve`](https://rdrr.io/r/stats/convolve.html)). The lag
  maximising cross-correlation is the optimal shift. Robust to
  multi-peak shapes.

- `"landmark"`:

  Like `"peak"`, but uses a user-supplied function
  `landmark.func(curve) -> index` to locate the feature of interest
  (e.g., steepest rise, first zero crossing).

- `"iterative"`:

  Alternate cross-correlation rotation and elastic alignment for up to
  `max.iter` iterations. Useful when the optimal rotation depends on the
  alignment and vice versa.

## References

Srivastava, A. and Klassen, E. (2016). *Functional and Shape Data
Analysis*. Springer.

## Examples

``` r
# Periodic sinusoidal curves with random circular shifts
argvals <- seq(0, 2 * pi, length.out = 100)
data <- matrix(0, 10, 100)
for (i in 1:10) {
  shift <- sample(0:99, 1)
  idx <- ((seq_len(100) - 1 + shift) %% 100) + 1
  data[i, ] <- sin(argvals[idx])
}
fd <- fdata(data, argvals = argvals)
rot <- periodic.rotate(fd)

# Cross-correlation method (robust to multi-peak curves)
rot_xcorr <- periodic.rotate(fd, method = "xcorr")
```
