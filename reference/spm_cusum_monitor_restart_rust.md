# CUSUM monitoring with restart after each alarm.

CUSUM monitoring with restart after each alarm.

## Usage

``` r
spm_cusum_monitor_restart_rust(
  sequential_data,
  argvals,
  rotation,
  mean,
  singular_values,
  centered,
  eigenvalues,
  t2_ucl,
  t2_alpha,
  t2_description,
  spe_ucl,
  spe_alpha,
  spe_description,
  chart_ncomp,
  config_alpha,
  config_tuning_fraction,
  config_seed,
  k,
  h,
  cusum_ncomp,
  cusum_alpha,
  multivariate
)
```
