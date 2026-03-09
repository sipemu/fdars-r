# Andrews Loadings: Project FPCA Eigenfunctions to Original Variables

Given FPCA results on Andrews-transformed data, projects each
eigenfunction back onto the Fourier basis to obtain loadings in the
original variable space. This reveals which variables contribute most to
each principal component.

## Usage

``` r
andrews_loadings(fpca, fd_andrews, ncomp = NULL)
```

## Arguments

- fpca:

  Result from
  [`fdata2pc()`](https://sipemu.github.io/fdars-r/reference/fdata2pc.md)
  (class `"fdata2pc"`).

- fd_andrews:

  The [fdata](https://sipemu.github.io/fdars-r/reference/fdata.md)
  object returned by
  [`andrews_transform()`](https://sipemu.github.io/fdars-r/reference/andrews_transform.md)
  (must carry the `"andrews_basis"` attribute).

- ncomp:

  Number of principal components to project (default: all available in
  `fpca`).

## Value

A data.frame of class `"andrews_loadings"` with columns:

- `Variable`:

  Original variable name.

- `Loading`:

  Inner-product loading value.

- `PC`:

  Principal component label (`"PC1"`, `"PC2"`, ...).

## See also

[`andrews_transform()`](https://sipemu.github.io/fdars-r/reference/andrews_transform.md),
[`fdata2pc()`](https://sipemu.github.io/fdars-r/reference/fdata2pc.md)

## Examples

``` r
X <- scale(iris[, 1:4])
fd <- andrews_transform(X)
fpca <- fdata2pc(fd, ncomp = 3)
loadings <- andrews_loadings(fpca, fd)
head(loadings)
#>       Variable     Loading  PC
#> 1 Sepal.Length -0.16355095 PC1
#> 2  Sepal.Width  0.08487338 PC1
#> 3 Petal.Length -0.18343483 PC1
#> 4  Petal.Width -0.17786725 PC1
#> 5 Sepal.Length  0.11957088 PC2
#> 6  Sepal.Width  0.29061451 PC2
```
