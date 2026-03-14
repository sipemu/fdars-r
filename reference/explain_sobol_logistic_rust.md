# Sobol indices for logistic model (no y, uses n_samples/seed)

Sobol indices for logistic model (no y, uses n_samples/seed)

## Usage

``` r
explain_sobol_logistic_rust(
  intercept,
  beta_t,
  beta_se,
  gamma,
  probabilities,
  predicted_classes,
  ncomp,
  accuracy,
  std_errors,
  coefficients,
  log_likelihood,
  iterations,
  fpca_mean,
  fpca_rotation_data,
  fpca_rotation_nrow,
  fpca_rotation_ncol,
  fpca_scores_data,
  fpca_scores_nrow,
  fpca_scores_ncol,
  data,
  n_samples,
  seed
)
```
