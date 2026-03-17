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
#>  [1] 0.6962850 0.1937026 0.5414913 0.2283806 0.6235426 0.5292162 0.1999323
#>  [8] 0.7119285 0.3644720 0.5518433 0.1967200 0.2944433 0.2669759 0.3795362
#> [15] 0.1396571 0.1349786 0.1595666 0.2356083 0.5955438 0.2597674 0.2412882
#> [22] 0.3711429 0.3422697 0.2319775 0.2597991 0.2445493 0.7257824 0.4593558
#> [29] 0.6573595 0.5711324 0.2698134 0.6638244 0.6730164 0.6028522 0.2569189
#> [36] 0.2872857 0.7505596 0.3923126 0.6604708 0.2058081 0.2391108 0.6252753
#> [43] 0.5794011 0.3457562 0.1684882 0.5865541 0.5920372 0.7775842 0.5776515
#> [50] 0.2548728
# }
```
