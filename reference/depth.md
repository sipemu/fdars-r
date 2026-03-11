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
#>  [1] 0.64 0.64 0.48 0.44 0.44 0.40 0.54 0.46 0.46 0.50
depth(fd, method = "mode")
#>  [1] 0.3111973 0.3146904 0.1598219 0.1540158 0.3041319 0.2406552 0.2954819
#>  [8] 0.2503295 0.2403409 0.2184980
depth(fd, method = "RP")
#>  [1] 0.3381818 0.3109091 0.2054545 0.1800000 0.3181818 0.2472727 0.3090909
#>  [8] 0.2836364 0.2600000 0.2745455
```
