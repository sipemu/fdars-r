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

## Examples

``` r
# \donttest{
set.seed(42)
t_grid <- seq(0, 1, length.out = 30)
X <- matrix(0, 40, 30)
for (i in 1:40) X[i, ] <- sin(2*pi*t_grid) * (2*(i > 20) - 1) + rnorm(30, sd = 0.2)
fd <- fdata(X, argvals = t_grid)
y_bin <- as.numeric(1:40 > 20)
fit <- functional.logistic(fd, y_bin)
result <- fregre.calibration(fit, as.numeric(as.character(y_bin)))
# }
```
