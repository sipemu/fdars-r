# TSRVF from Pre-computed Alignment

Compute the TSRVF representation from a pre-computed Karcher mean,
avoiding the expensive re-computation of the alignment.

## Usage

``` r
tsrvf.from.alignment(karcher)
```

## Arguments

- karcher:

  An object of class 'karcher.mean'.

## Value

An object of class 'tsrvf'.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
km <- karcher.mean(fd, max.iter = 5)
tv <- tsrvf.from.alignment(km)
# }
```
