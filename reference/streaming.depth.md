# Streaming Depth for Functional Data

Compute depth values using a streaming (pre-sorted reference) algorithm.
This approach pre-sorts the reference data at each time point, enabling
O(T log N) depth computation per curve instead of O(N^2 T).

## Usage

``` r
streaming.depth(fdataobj, fdataori = NULL, method = c("FM", "MBD", "BD"), ...)
```

## Arguments

- fdataobj:

  An object of class 'fdata' to compute depth for.

- fdataori:

  An object of class 'fdata' as reference sample. If NULL, uses fdataobj
  as reference.

- method:

  Streaming depth method: "FM" (Fraiman-Muniz, default), "MBD" (modified
  band depth), or "BD" (band depth).

- ...:

  Additional arguments (ignored).

## Value

Numeric vector of depth values, one per curve in fdataobj.

## References

Fraiman, R. and Muniz, G. (2001). Trimmed means for functional data.
*Test*, 10(2):419–440.

Lopez-Pintado, S. and Romo, J. (2009). On the concept of depth for
functional data. *Journal of the American Statistical Association*,
104(486):718–734.

## Examples

``` r
fd <- fdata(matrix(rnorm(200), 20, 10))
streaming.depth(fd, method = "FM")
#>  [1] 0.41 0.37 0.45 0.49 0.62 0.54 0.56 0.59 0.50 0.45 0.53 0.54 0.49 0.39 0.52
#> [16] 0.57 0.60 0.37 0.52 0.49
```
