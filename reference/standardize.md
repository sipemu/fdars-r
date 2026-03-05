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
#>  [1]  4.440892e-17 -2.692291e-16  6.591949e-18 -5.134781e-17 -3.580469e-16
#>  [6] -2.081668e-18  2.192690e-16  2.026157e-16  3.415671e-16 -7.494005e-17
apply(fd_std$data, 1, sd)
#>  [1] 1 1 1 1 1 1 1 1 1 1
```
