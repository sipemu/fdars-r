# Depth Functions for Functional Data

Functions for computing various depth measures for functional data.
Compute Functional Data Depth

## Usage

``` r
depth(
  fdataobj,
  fdataori = NULL,
  method = c("FM", "mode", "RP", "RT", "BD", "MBD", "MEI", "FSD", "KFSD", "RPD",
    "streaming"),
  ...
)
```

## Arguments

- fdataobj:

  An object of class 'fdata' to compute depth for.

- fdataori:

  An object of class 'fdata' as reference sample. If NULL, uses fdataobj
  as reference.

- method:

  Depth method to use. One of "FM" (Fraiman-Muniz), "mode" (modal), "RP"
  (random projection), "RT" (random Tukey), "BD" (band depth), "MBD"
  (modified band depth), "MEI" (modified epigraph index), "FSD"
  (functional spatial depth), "KFSD" (kernel functional spatial depth),
  "RPD" (random projection with derivatives), or "streaming" (streaming
  depth via pre-sorted reference). Default is "FM".

- ...:

  Additional arguments passed to the specific depth function.

## Value

A numeric vector of depth values, one per curve in fdataobj.

## Details

Unified interface for computing various depth measures for functional
data.

Available methods:

- FM:

  Fraiman-Muniz depth - integrates univariate depths over domain

- mode:

  Modal depth - based on kernel density estimation

- RP:

  Random projection depth - projects to random directions

- RT:

  Random Tukey depth - halfspace depth via random projections

- BD:

  Band depth - proportion of bands containing the curve (1D only)

- MBD:

  Modified band depth - allows partial containment (1D only)

- MEI:

  Modified epigraph index - proportion of time below other curves (1D
  only)

- FSD:

  Functional spatial depth - based on spatial signs

- KFSD:

  Kernel functional spatial depth - smoothed FSD

- RPD:

  Random projection with derivatives - includes curve derivatives

- streaming:

  Streaming depth via pre-sorted reference (FM/MBD/BD selectable via
  `streaming.method`)

## Examples

``` r
fd <- fdata(matrix(rnorm(100), 10, 10))

# Different depth methods
depth(fd, method = "FM")
#>  [1] 0.60 0.46 0.42 0.60 0.46 0.44 0.50 0.44 0.50 0.58
depth(fd, method = "mode")
#>  [1] 0.2631642 0.2007614 0.1992753 0.2377022 0.1550116 0.1902275 0.2261796
#>  [8] 0.1947884 0.2284591 0.2331016
depth(fd, method = "RP")
#>  [1] 0.3218182 0.2872727 0.2509091 0.2872727 0.2418182 0.2418182 0.2872727
#>  [8] 0.2454545 0.2563636 0.3072727
```
