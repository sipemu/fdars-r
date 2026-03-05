# Compute GCV score for basis fit GCV = RSS/n / (1 - edf/n)^2 When pooled=true: compute single GCV across all curves When pooled=false: compute per-curve GCV and return mean

Compute GCV score for basis fit GCV = RSS/n / (1 - edf/n)^2 When
pooled=true: compute single GCV across all curves When pooled=false:
compute per-curve GCV and return mean

## Usage

``` r
basis_gcv_1d(data, argvals, nbasis, basis_type, lambda, pooled)
```
