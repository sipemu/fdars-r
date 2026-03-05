# Detect amplitude modulation in seasonal time series using Hilbert transform

Detect amplitude modulation in seasonal time series using Hilbert
transform

## Usage

``` r
seasonal_detect_amplitude_modulation(
  data,
  argvals,
  period,
  modulation_threshold,
  seasonality_threshold
)
```

## Arguments

- data:

  Matrix of functional data (n x m)

- argvals:

  Vector of evaluation points (length m)

- period:

  Seasonal period in argvals units

- modulation_threshold:

  CV threshold for detecting modulation (default: 0.15)

- seasonality_threshold:

  Strength threshold for seasonality (default: 0.3)

## Value

List with detection results
