# L1 (Least Absolute Deviation) Functional Regression

Fits a functional linear model using L1 loss (median regression), which
is robust to outliers in the response.

## Usage

``` r
fregre.l1(fdataobj, y, scalar.covariates = NULL, ncomp = 3)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- y:

  Response vector.

- scalar.covariates:

  Optional matrix of scalar covariates.

- ncomp:

  Number of FPC components (default 3).

## Value

An object of class 'fregre.robust' with components:

- intercept:

  Intercept

- beta.t:

  Functional coefficient beta(t)

- fitted.values:

  Fitted response values

- residuals:

  Residuals

- coefficients:

  Regression coefficients

- ncomp:

  Number of FPC components used

- iterations:

  Number of IRLS iterations

- converged:

  Whether the algorithm converged

- weights:

  Final IRLS weights

- r.squared:

  R-squared statistic

- method:

  Either "l1" or "huber"

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.l1(fd, y, ncomp = 3)
# }
```
