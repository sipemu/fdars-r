# Penalized Basis Smoothing with GCV-Optimal Lambda

Selects the smoothing parameter by minimizing the generalized
cross-validation (GCV) criterion over a grid of lambda values.

## Usage

``` r
smooth.basis.gcv(
  fdataobj,
  type = c("bspline", "fourier"),
  nbasis = NULL,
  lfd.order = 2,
  log.lambda.range = c(-10, 5),
  n.grid = 50,
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

- lfd.order:

  Order of the linear differential operator for the penalty.

- log.lambda.range:

  Range of log10(lambda) to search over.

- n.grid:

  Number of grid points for the lambda search.

- period:

  Period for Fourier basis.

## Value

An object of class 'smooth.basis' (same as smooth.basis.fd).

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
sm <- smooth.basis.gcv(fd, type = "bspline", nbasis = 5)
sm
#> Penalized Basis Smoothing
#>   Basis type: bspline 
#>   Number of basis functions: 6 
#>   Lambda: 15.64 
#>   EDF: 2.1 
#>   GCV: 1.19421 
# }
```
