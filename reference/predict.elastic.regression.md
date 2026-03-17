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
#>  [1] -3.11306228 -1.42688868 -1.61107456  0.74311441 -1.06705498  1.18541705
#>  [7] -2.09651172 -2.56909387  1.84083574 -1.11636204  1.63017188  1.21584324
#> [13] -0.97711056  0.35245738  3.65096690 -0.58910284 -2.13937372  1.73379774
#> [19] -0.54451970 -2.61876216 -1.55627304 -2.31276349 -2.16419619 -0.04063609
#> [25]  0.42662450 -1.40038505  0.24870179  1.67840446 -0.61629839  0.56299841
#> [31]  0.86709024  1.05335288  1.58092182  2.86459995  2.52631782  1.47876936
#> [37] -1.73796221 -0.68022813 -1.41161398 -0.56503697  0.84710761 -0.41294739
#> [43] -0.74856594  0.57063195  0.54395287 -1.47500367 -1.35572821  0.91602643
#> [49] -1.78592239 -1.33549528
# }
```
