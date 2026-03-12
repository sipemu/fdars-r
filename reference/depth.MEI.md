# Modified Epigraph Index

Wrapper for `depth(method = "MEI")`. Useful as a function argument.

## Usage

``` r
depth.MEI(fdataobj, fdataori = NULL, ...)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- fdataori:

  Reference sample (default: fdataobj itself).

- ...:

  Additional arguments.

## Value

Numeric vector of depth values.

## See also

[`depth`](https://sipemu.github.io/fdars-r/reference/depth.md)

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(50), nrow = 5), argvals = seq(0, 1, length.out = 10))
d <- depth.MEI(fd)
# }
```
