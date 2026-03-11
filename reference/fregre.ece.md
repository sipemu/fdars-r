# Expected Calibration Error (Logistic)

Computes ECE, MCE, and ACE for a functional logistic model.

## Usage

``` r
fregre.ece(model, y, n.bins = 10)
```

## Arguments

- model:

  A fitted functional logistic model.

- y:

  True binary labels (0/1).

- n.bins:

  Number of bins (default 10).

## Value

A list with ece, mce, ace, n_bins, and bin_ece_contributions.
