# High-Accuracy Gradient for Functional Data

Computes the numerical gradient (first derivative) of functional data
using a 5-point stencil for uniform grids (O(h^4) accuracy) or 3-point
Lagrange for non-uniform grids. Complementary to `deriv.fdata` which
uses central differences.

## Usage

``` r
fdata.gradient(fdataobj)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

## Value

An `fdata` object with the gradient values.

## Examples

``` r
t <- seq(0, 1, length.out = 100)
X <- matrix(sin(2 * pi * t), nrow = 1)
fd <- fdata(X, argvals = t)
gd <- fdata.gradient(fd)
```
