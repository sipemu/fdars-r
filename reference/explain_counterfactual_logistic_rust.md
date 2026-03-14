# Counterfactual for logistic model (max_iter/step_size instead of target_value)

Counterfactual for logistic model (max_iter/step_size instead of
target_value)

## Usage

``` r
explain_counterfactual_logistic_rust(
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
  observation,
  max_iter,
  step_size
)
```
