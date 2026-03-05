# Compute AIC for basis fit AIC = n \* log(RSS/n) + 2 \* total_edf Where total_edf = n_curves \* edf (each curve has edf parameters) When pooled=true: compute single AIC across all curves When pooled=false: compute per-curve AIC and return mean

Compute AIC for basis fit AIC = n \* log(RSS/n) + 2 \* total_edf Where
total_edf = n_curves \* edf (each curve has edf parameters) When
pooled=true: compute single AIC across all curves When pooled=false:
compute per-curve AIC and return mean

## Usage

``` r
basis_aic_1d(data, argvals, nbasis, basis_type, lambda, pooled)
```
