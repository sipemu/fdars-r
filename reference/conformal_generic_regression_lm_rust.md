# Generic conformal regression for a pre-fitted fregre.lm model

Generic conformal regression for a pre-fitted fregre.lm model

## Usage

``` r
conformal_generic_regression_lm_rust(
  data,
  y,
  test_data,
  scalar_train,
  scalar_test,
  fpca_mean,
  fpca_rotation,
  fpca_scores,
  ncomp,
  coefficients,
  gamma,
  calibration_indices,
  cal_fraction,
  alpha,
  seed
)
```
