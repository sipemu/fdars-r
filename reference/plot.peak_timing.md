# Plot method for peak_timing objects

Plots peak timing (day within cycle) and peak value across cycles, with
a linear trend line showing whether peaks shift earlier or later.

## Usage

``` r
# S3 method for class 'peak_timing'
plot(x, period = NULL, ...)
```

## Arguments

- x:

  A peak_timing object (from
  [`analyze.peak.timing`](https://sipemu.github.io/fdars-r/reference/analyze.peak.timing.md)).

- period:

  The period used for analysis (for computing day-of-cycle).

- ...:

  Additional arguments (ignored).

## Value

A ggplot object.
