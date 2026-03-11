# Significant Regions

Identifies domain regions where the beta coefficient is significantly
different from zero based on standard errors.

## Usage

``` r
fregre.significant.regions(model, alpha = 0.05)
```

## Arguments

- model:

  A fitted `fregre.lm` model.

- alpha:

  Significance level (default 0.05).

## Value

A list with start_idx, end_idx, and direction vectors, or NULL if no
significant regions found.
