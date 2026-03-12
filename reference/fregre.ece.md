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

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y_bin <- factor(sample(0:1, 50, replace = TRUE))
fit <- functional.logistic(fd, y_bin)
#> Error in functional.logistic(fd, y_bin): functional.logistic failed: check data dimensions and response
result <- fregre.ece(fit, as.numeric(as.character(y_bin)))
#> Error: object 'fit' not found
# }
```
