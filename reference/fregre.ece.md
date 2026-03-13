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
set.seed(42)
t_grid <- seq(0, 1, length.out = 30)
X <- matrix(0, 40, 30)
for (i in 1:40) X[i, ] <- sin(2*pi*t_grid) * (2*(i > 20) - 1) + rnorm(30, sd = 0.2)
fd <- fdata(X, argvals = t_grid)
y_bin <- as.numeric(1:40 > 20)
fit <- functional.logistic(fd, y_bin)
result <- fregre.ece(fit, as.numeric(as.character(y_bin)))
# }
```
