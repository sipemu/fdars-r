# Adaptive EWMA (AMEWMA) monitoring on sequential functional data.

Adaptive EWMA (AMEWMA) monitoring on sequential functional data.

## Usage

``` r
spm_amewma_monitor_rust(
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
  lambda_min,
  lambda_max,
  lambda_init,
  eta,
  amewma_ncomp,
  amewma_alpha
)
```
