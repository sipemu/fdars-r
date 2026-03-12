# Streaming Depth (Alias)

Alias for
[`streaming.depth`](https://sipemu.github.io/fdars-r/reference/streaming.depth.md).

## Usage

``` r
depth.streaming(fdataobj, fdataori = NULL, ...)
```

## Arguments

- fdataobj:

  An object of class 'fdata' to compute depth for.

- fdataori:

  An object of class 'fdata' as reference sample. If NULL, uses fdataobj
  as reference.

- ...:

  Additional arguments (ignored).

## Value

Numeric vector of depth values.

## See also

[`streaming.depth`](https://sipemu.github.io/fdars-r/reference/streaming.depth.md),
[`depth`](https://sipemu.github.io/fdars-r/reference/depth.md)

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(50), nrow = 5), argvals = seq(0, 1, length.out = 10))
d <- depth.streaming(fd)
# }
```
