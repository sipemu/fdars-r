# Compute BIC for basis fit BIC = n \* log(RSS/n) + log(n) \* total_edf Where total_edf = n_curves \* edf (each curve has edf parameters) When pooled=true: compute single BIC across all curves When pooled=false: compute per-curve BIC and return mean

Compute BIC for basis fit BIC = n \* log(RSS/n) + log(n) \* total_edf
Where total_edf = n_curves \* edf (each curve has edf parameters) When
pooled=true: compute single BIC across all curves When pooled=false:
compute per-curve BIC and return mean

## Usage

``` r
basis_bic_1d(data, argvals, nbasis, basis_type, lambda, pooled)
```
