# Periodic Rotation for Functional Data

Circularly rotate each curve so that its global maximum is at a
canonical grid position. This is useful as a pre-processing step before
elastic alignment of periodic functional data (e.g., data on \\\[0,
2\pi\]\\ where \\f(0) = f(2\pi)\\).

## Usage

``` r
periodic.rotate(fdataobj, target.pos = NULL)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- target.pos:

  Target grid index for the global maximum. If NULL (default), uses
  `floor(ncol(fdataobj$data) / 4)`.

## Value

A list with components:

- fdataobj:

  fdata of rotated curves

- shifts:

  integer vector of circular shifts applied to each curve

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
```
