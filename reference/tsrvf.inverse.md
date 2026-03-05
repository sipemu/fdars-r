# Inverse TSRVF Transform

Reconstruct aligned curves from their TSRVF tangent vector
representation.

## Usage

``` r
tsrvf.inverse(tsrvf.obj)
```

## Arguments

- tsrvf.obj:

  An object of class 'tsrvf'.

## Value

An fdata object of reconstructed aligned curves.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
tv <- tsrvf.transform(fd, max.iter = 5)
recon <- tsrvf.inverse(tv)
# }
```
