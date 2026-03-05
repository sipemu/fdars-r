# Landmark Registration for Functional Data

Feature-based alignment using detected landmarks (peaks, valleys,
zero-crossings, inflection points). Detect Landmarks in Functional Data

## Usage

``` r
detect.landmarks(
  fdataobj,
  kind = c("peak", "valley", "zero", "inflection"),
  min.prominence = 0
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- kind:

  Character vector of landmark types to detect. One of "peak", "valley",
  "zero", "inflection". Default is "peak".

- min.prominence:

  Minimum prominence for peak/valley detection. Default is 0.

## Value

A list of data frames (one per curve), each with columns: position,
kind, value, prominence.

## Details

Detect feature landmarks (peaks, valleys, zero-crossings, inflection
points) in each curve of a functional data object.

## Examples

``` r
t <- seq(0, 1, length.out = 100)
X <- matrix(0, 5, 100)
for (i in 1:5) X[i, ] <- sin(2*pi*t + i/10)
fd <- fdata(X, argvals = t)
lms <- detect.landmarks(fd, kind = "peak")
```
