# Predict from Elastic Logistic Classification

Predict from Elastic Logistic Classification

## Usage

``` r
# S3 method for class 'elastic.logistic'
predict(object, newdata, type = c("response", "class"), ...)
```

## Arguments

- object:

  A fitted object of class 'elastic.logistic'.

- newdata:

  New functional data (fdata object).

- type:

  "response" for probabilities (default) or "class" for classes.

- ...:

  Additional arguments (ignored).

## Value

Predicted probabilities or class labels.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- sample(0:1, 50, replace = TRUE)
fit <- elastic.logistic(fd, y)
predict(fit, fd)
#>  [1] 0.23395122 0.58746347 0.25313689 0.09666869 0.29484102 0.27019420
#>  [7] 0.11468709 0.26701090 0.57734466 0.19549513 0.21835777 0.25696519
#> [13] 0.26092479 0.31485759 0.24869446 0.79157420 0.18371640 0.56067194
#> [19] 0.25803489 0.46109036 0.08593395 0.61388580 0.47692456 0.20833996
#> [25] 0.21849399 0.24659935 0.25419470 0.31282411 0.14228280 0.13737015
#> [31] 0.57698349 0.42497694 0.32369637 0.48529685 0.27738131 0.09587824
#> [37] 0.52723364 0.14259232 0.19467879 0.54051360 0.36949023 0.15658381
#> [43] 0.17856707 0.53765959 0.45595977 0.47049241 0.57381204 0.21374388
#> [49] 0.30597499 0.18382864
# }
```
