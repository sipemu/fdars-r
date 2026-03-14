# Conditional permutation importance for logistic model

Conditional permutation importance for logistic model

## Usage

``` r
explain_conditional_importance_logistic_rust(
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
  y,
  n_bins,
  n_perm,
  seed
)
```
