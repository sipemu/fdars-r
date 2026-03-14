# Huber M-Estimation Functional Regression

Fits a functional linear model using Huber's M-estimation, which
provides a smooth compromise between least squares and L1.

## Usage

``` r
fregre.huber(fdataobj, y, scalar.covariates = NULL, ncomp = 3, k = 1.345)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- y:

  Response vector.

- scalar.covariates:

  Optional matrix of scalar covariates.

- ncomp:

  Number of FPC components (default 3).

- k:

  Huber tuning parameter (default 1.345 for 95 percent efficiency).

## Value

An object of class 'fregre.robust'. See
[`fregre.l1`](https://sipemu.github.io/fdars-r/reference/fregre.l1.md).

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.huber(fd, y, ncomp = 3)
# }
```
