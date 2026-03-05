# Elastic Phase-Amplitude Decomposition

Decompose the difference between two curves into amplitude and phase
components using the elastic (Fisher-Rao) framework.

## Usage

``` r
elastic.decomposition(f1, f2, argvals = NULL, lambda = 0)
```

## Arguments

- f1:

  An 'fdata' object with a single curve, or a numeric vector.

- f2:

  An 'fdata' object with a single curve, or a numeric vector.

- argvals:

  Evaluation points (required if f1/f2 are numeric vectors).

- lambda:

  Penalty weight on warp deviation from identity. Default 0.

## Value

A list with components:

- gamma:

  Warping function

- f_aligned:

  Aligned version of f2

- d_amplitude:

  Amplitude distance

- d_phase:

  Phase distance

## Examples

``` r
# \donttest{
t <- seq(0, 1, length.out = 50)
f1 <- fdata(matrix(sin(2*pi*t), 1, 50), argvals = t)
f2 <- fdata(matrix(sin(2*pi*t + 0.3), 1, 50), argvals = t)
dec <- elastic.decomposition(f1, f2)
# }
```
