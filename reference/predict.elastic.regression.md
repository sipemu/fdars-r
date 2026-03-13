# Predict from Elastic Regression

Predict from Elastic Regression

## Usage

``` r
# S3 method for class 'elastic.regression'
predict(object, newdata, ...)
```

## Arguments

- object:

  A fitted object of class 'elastic.regression'.

- newdata:

  New functional data (fdata object).

- ...:

  Additional arguments (ignored).

## Value

Predicted response values.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- elastic.regression(fd, y)
predict(fit, fd)
#>  [1] -1.00094012  2.23323783 -1.78061770 -1.01138748  3.47114563  2.13518956
#>  [7] -1.93812711  1.13792798  1.62684810  0.64315401  1.97236954 -0.58515802
#> [13]  1.06372947  0.70535918  0.09716373 -2.07121537 -0.27912513 -0.86858471
#> [19] -1.18767692  0.58454993 -1.64405677 -0.11711361  0.31386114  1.59557050
#> [25] -0.86971525  0.97318034  1.19988522 -1.26552440 -0.11636136 -1.20088597
#> [31]  1.49071253 -1.27958344  0.41825813  1.90561462  0.74744636  0.49060561
#> [37] -2.20762327 -0.44469409  2.15683069  0.15793466  1.38909360 -2.97324543
#> [43]  0.73990995  0.61585453 -0.98159381  0.64654942 -1.37299511  0.50244021
#> [49]  1.09736047  1.41740013
# }
```
