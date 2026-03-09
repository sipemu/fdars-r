# Andrews Transformation

Transforms a multivariate dataset into functional data using the Andrews
Fourier expansion. Each row of the input matrix becomes a curve on
\\\[-\pi, \pi\]\\ defined by \$\$f_i(t) = x\_{i1}/\sqrt{2} +
x\_{i2}\sin(t) + x\_{i3}\cos(t) + x\_{i4}\sin(2t) + x\_{i5}\cos(2t) +
\ldots\$\$

## Usage

``` r
andrews_transform(X, m = 200)
```

## Arguments

- X:

  Numeric matrix or data.frame (\\n \times p\\). Rows are observations,
  columns are variables. Coerced via
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html).

- m:

  Number of grid points on \\\[-\pi, \pi\]\\ (default 200).

## Value

An [fdata](https://sipemu.github.io/fdars-r/reference/fdata.md) object
(\\n\\ curves, \\m\\ grid points) with additional attributes:

- `andrews_basis`:

  The \\m \times p\\ Fourier basis matrix.

- `andrews_varnames`:

  Variable names (`colnames(X)` or `paste0("V", 1:p)` if `NULL`).

## Details

The resulting `fdata` object carries the Fourier basis matrix and
variable names as attributes, enabling
[`andrews_loadings()`](https://sipemu.github.io/fdars-r/reference/andrews_loadings.md)
to project FPCA eigenfunctions back to the original variable space.

## See also

[`andrews_loadings()`](https://sipemu.github.io/fdars-r/reference/andrews_loadings.md),
[`fdata2pc()`](https://sipemu.github.io/fdars-r/reference/fdata2pc.md)

## Examples

``` r
X <- scale(iris[, 1:4])
fd <- andrews_transform(X)
plot(fd)
```
