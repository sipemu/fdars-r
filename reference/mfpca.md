# Multivariate Functional Principal Component Analysis

Extends univariate FPCA to handle multiple functional variables observed
on potentially different grids. Variables are optionally weighted by
their inverse standard deviation before joint SVD.

## Usage

``` r
mfpca(fdataobj.list, ncomp = 5, weighted = TRUE)
```

## Arguments

- fdataobj.list:

  A list of `fdata` objects, one per functional variable. All must have
  the same number of observations (rows) but may differ in grid size.

- ncomp:

  Number of principal components to extract (default 5).

- weighted:

  Logical; whether to weight each variable by 1/std_dev before SVD
  (default TRUE).

## Value

An object of class `mfpca` with components:

- scores:

  Score matrix (n x ncomp)

- eigenfunctions:

  List of eigenfunction matrices, one per variable

- eigenvalues:

  Eigenvalues (length ncomp)

- means:

  List of mean functions, one per variable

- scales:

  Per-variable standard deviations

- grid.sizes:

  Grid sizes per variable

## Examples

``` r
# \donttest{
# Two functional variables
n <- 40
fd1 <- fdata(matrix(rnorm(n * 20), n, 20), argvals = seq(0, 1, length.out = 20))
fd2 <- fdata(matrix(rnorm(n * 15), n, 15), argvals = seq(0, 1, length.out = 15))
result <- mfpca(list(fd1, fd2), ncomp = 3)
dim(result$scores)
#> [1] 40  3
# }
```
