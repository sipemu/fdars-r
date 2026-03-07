# FPC-based Function-on-Scalar Regression

Fits a function-on-scalar regression using functional principal
components.

## Usage

``` r
fosr.fpc(fdataobj, predictors, ncomp = 3)
```

## Arguments

- fdataobj:

  An object of class 'fdata' (functional response).

- predictors:

  A matrix of scalar predictors (n x p).

- ncomp:

  Number of FPC components (default 3).

## Value

An object of class 'fosr'.
