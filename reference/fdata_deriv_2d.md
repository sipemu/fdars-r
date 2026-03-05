# Compute 2D partial derivatives for surface data

For a surface f(s,t), computes:

- ds: partial derivative with respect to s (∂f/∂s)

- dt: partial derivative with respect to t (∂f/∂t)

- dsdt: mixed partial derivative (∂²f/∂s∂t)

## Usage

``` r
fdata_deriv_2d(data, argvals_s, argvals_t, m1, m2)
```

## Details

Data layout: n surfaces, each stored as m1*m2 values in row-major order
(s varies fastest) Return: list with three matrices (ds, dt, dsdt), each
n x (m1*m2)
