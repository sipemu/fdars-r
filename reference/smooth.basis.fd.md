# Penalized Basis Smoothing

Smooth functional data using penalized B-spline or Fourier basis
expansion. Penalized Basis Smoothing

## Usage

``` r
smooth.basis.fd(
  fdataobj,
  type = c("bspline", "fourier"),
  nbasis = NULL,
  lambda = 0,
  lfd.order = 2,
  period = NULL
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- type:

  Basis type: "bspline" or "fourier".

- nbasis:

  Number of basis functions. If NULL, defaults to min(ncol/4, 30).

- lambda:

  Smoothing parameter (penalty weight). Use 0 for no penalty.

- lfd.order:

  Order of the linear differential operator for the penalty (default 2 =
  penalize curvature).

- period:

  Period for Fourier basis (only used when type = "fourier"). If NULL,
  defaults to the range of argvals.

## Value

An object of class 'smooth.basis' with components:

- coefficients:

  Matrix of basis coefficients (nbasis x n)

- fitted:

  Fitted smoothed curves as an fdata object

- edf:

  Effective degrees of freedom

- gcv:

  Generalized cross-validation score

- aic:

  Akaike information criterion

- bic:

  Bayesian information criterion

- nbasis:

  Number of basis functions used

- lambda:

  Smoothing parameter used

- type:

  Basis type used

- fdataobj:

  Original functional data

## Details

Fits a penalized basis expansion to smooth functional data. Supports
B-spline and Fourier basis types with a roughness penalty.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
sm <- smooth.basis.fd(fd, type = "bspline", nbasis = 5, lambda = 0.1)
sm
#> Penalized Basis Smoothing
#>   Basis type: bspline 
#>   Number of basis functions: 6 
#>   Lambda: 0.1 
#>   EDF: 2.2 
#>   GCV: 1.31118 
# }
```
