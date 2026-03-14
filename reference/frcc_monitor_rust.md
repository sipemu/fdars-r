# FRCC Phase II: Monitor new data against a functional regression control chart.

FRCC Phase II: Monitor new data against a functional regression control
chart.

## Usage

``` r
frcc_monitor_rust(
  new_y,
  new_predictors,
  argvals,
  fosr_intercept,
  fosr_beta,
  fosr_lambda,
  resid_rotation,
  resid_mean,
  resid_singular_values,
  resid_centered,
  eigenvalues,
  t2_ucl,
  t2_alpha,
  t2_description,
  spe_ucl,
  spe_alpha,
  spe_description,
  ncomp,
  config_fosr_lambda,
  config_alpha,
  config_tuning_fraction,
  config_seed
)
```
