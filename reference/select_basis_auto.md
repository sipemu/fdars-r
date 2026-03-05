# Automatic basis selection for each curve individually.

This function compares Fourier and P-spline bases for each curve,
selecting the optimal basis type and number of basis functions using
model selection criteria (GCV, AIC, or BIC).

## Usage

``` r
select_basis_auto(
  data,
  argvals,
  criterion,
  nbasis_min,
  nbasis_max,
  lambda_pspline,
  use_seasonal_hint
)
```

## Arguments

- `data` - Functional data matrix (n x m)

- `argvals` - Evaluation points

- `criterion` - Model selection criterion: 0=GCV (default), 1=AIC, 2=BIC

- `nbasis_min` - Minimum number of basis functions (0 for auto)

- `nbasis_max` - Maximum number of basis functions (0 for auto)

- `lambda_pspline` - Smoothing parameter for P-spline (negative for
  auto-select)

- `use_seasonal_hint` - Whether to use FFT to detect seasonality

## Returns

Named list with:

- basis_type: Integer vector (0=pspline, 1=fourier)

- nbasis: Integer vector of selected nbasis per curve

- score: Numeric vector of criterion scores

- coefficients: List of coefficient vectors

- fitted: Fitted values matrix (n x m)

- edf: Numeric vector of effective degrees of freedom

- seasonal_detected: Logical vector

- lambda: Numeric vector of lambda values (NA for Fourier)

- criterion: Character (criterion name used)
