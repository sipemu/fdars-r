# Landmark Registration

Register (align) functional data by detecting landmarks and warping
curves so that the landmarks are mapped to common target positions.

## Usage

``` r
landmark.register(
  fdataobj,
  kind = c("peak", "valley", "zero", "inflection"),
  min.prominence = 0,
  expected.count = 0
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- kind:

  Character: landmark type to detect. One of "peak", "valley", "zero",
  "inflection". Default is "peak".

- min.prominence:

  Minimum prominence for peak/valley detection. Default is 0.

- expected.count:

  Expected number of landmarks per curve. If 0 (default), uses the
  minimum number detected across all curves.

## Value

An object of class 'landmark.register' with components:

- registered:

  fdata of registered curves

- gammas:

  fdata of warping functions

- landmarks:

  list of detected landmarks per curve

- target_landmarks:

  common target landmark positions

- fdataobj:

  the original input

## Examples

``` r
# \donttest{
t <- seq(0, 1, length.out = 100)
X <- matrix(0, 10, 100)
for (i in 1:10) X[i, ] <- sin(2*pi*(t - i/50))
fd <- fdata(X, argvals = t)
lr <- landmark.register(fd, kind = "peak", min.prominence = 0.5)
# }
```
