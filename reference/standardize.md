# Standardize functional data (z-score normalization)

Transforms each curve to have mean 0 and standard deviation 1. This is
useful for comparing curve shapes regardless of their level or scale.

## Usage

``` r
standardize(fdataobj)

# S3 method for class 'fdata'
standardize(fdataobj)

# S3 method for class 'irregFdata'
standardize(fdataobj)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

## Value

A standardized 'fdata' object where each curve has mean 0 and sd 1.

## Examples

``` r
fd <- fdata(matrix(rnorm(100) * 10 + 50, 10, 10), argvals = seq(0, 1, length.out = 10))
fd_std <- standardize(fd)
# Check: each curve now has mean ~0 and sd ~1
rowMeans(fd_std$data)
#>  [1]  4.218847e-16 -3.552714e-16 -9.992007e-17 -8.881784e-17 -1.110223e-15
#>  [6]  2.831069e-16 -5.107026e-16  1.065814e-15 -8.881784e-17 -5.107026e-16
apply(fd_std$data, 1, sd)
#>  [1] 1 1 1 1 1 1 1 1 1 1
```
