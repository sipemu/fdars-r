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
#>  [1] 0.56 0.42 0.52 0.50 0.62 0.48 0.32 0.62 0.62 0.34
depth(fd, method = "mode")
#>  [1] 0.2770756 0.1942053 0.2895437 0.3010254 0.3157327 0.2176444 0.1454690
#>  [8] 0.2934605 0.1814225 0.2107882
depth(fd, method = "RP")
#>  [1] 0.3072727 0.2145455 0.3436364 0.3127273 0.2727273 0.2745455 0.2145455
#>  [8] 0.3000000 0.2400000 0.2472727
```
