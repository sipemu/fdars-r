# Select Number of Principal Components

Selects the number of principal components to retain based on
eigenvalues using one of several criteria: cumulative variance
threshold, Kaiser rule, elbow method, or a fixed number.

## Usage

``` r
spm.ncomp.select(
  eigenvalues,
  method = c("variance90", "kaiser", "elbow", "fixed"),
  threshold = 0.9
)
```

## Arguments

- eigenvalues:

  Numeric vector of eigenvalues from FPCA.

- method:

  Selection method. One of `"variance90"` (cumulative variance exceeds
  `threshold`), `"kaiser"` (eigenvalues \> mean), `"elbow"` (scree plot
  elbow), or `"fixed"` (use `threshold` as the number of components).

- threshold:

  For `"variance90"`, the cumulative proportion of variance to exceed
  (default 0.9). For `"fixed"`, the exact number of components.

## Value

An integer: the selected number of components.

## Examples

``` r
# \donttest{
eigenvalues <- c(5.0, 2.0, 1.0, 0.5, 0.2, 0.1)
spm.ncomp.select(eigenvalues, method = "variance90", threshold = 0.9)
#> [1] 3
spm.ncomp.select(eigenvalues, method = "kaiser")
#> [1] 2
# }
```
