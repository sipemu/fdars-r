# Measure seasonal strength using wavelet (Morlet) method

Uses Continuous Wavelet Transform with Morlet wavelet to measure power
at the specified seasonal period.

## Usage

``` r
seasonal_strength_wavelet(data, argvals, period)
```

## Arguments

- data:

  Matrix of functional data (n x m)

- argvals:

  Vector of evaluation points (length m)

- period:

  Seasonal period in argvals units

## Value

Seasonal strength as ratio of wavelet power to total variance (0 to 1)
