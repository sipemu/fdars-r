# Calibration Diagnostics (Logistic)

Computes calibration diagnostics for a functional logistic model: Brier
score, log loss, and Hosmer-Lemeshow test.

## Usage

``` r
fregre.calibration(model, y, n.groups = 10)
```

## Arguments

- model:

  A fitted functional logistic model.

- y:

  True binary labels (0/1).

- n.groups:

  Number of groups for reliability diagram (default 10).

## Value

A list with brier_score, log_loss, hosmer_lemeshow_chi2,
hosmer_lemeshow_df, n_groups, reliability_bins, and bin_counts.
