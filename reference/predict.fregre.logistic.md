# Predict from Functional Logistic Model

Predict from Functional Logistic Model

## Usage

``` r
# S3 method for class 'fregre.logistic'
predict(
  object,
  newdata = NULL,
  new.scalar = NULL,
  type = c("response", "class"),
  ...
)
```

## Arguments

- object:

  A fitted object of class 'fregre.logistic'.

- newdata:

  New functional data (fdata object). If NULL, returns fitted
  probabilities.

- new.scalar:

  Optional matrix of new scalar covariates.

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
y <- rbinom(50, 1, 0.5)
fit <- functional.logistic(fd, y, ncomp = 3)
predict(fit)
#>  [1] 0.2157185 0.4495960 0.1840527 0.4529418 0.4040873 0.4404098 0.3415921
#>  [8] 0.4850081 0.4563026 0.3113535 0.5408615 0.3619769 0.4199375 0.7143978
#> [15] 0.2973238 0.4702335 0.3796918 0.5149536 0.5771892 0.4025855 0.2758588
#> [22] 0.6852547 0.4485272 0.6462660 0.5125092 0.3986676 0.4992176 0.4000180
#> [29] 0.2988760 0.5822917 0.4079047 0.3665283 0.4576293 0.4109688 0.5532176
#> [36] 0.4257777 0.1918132 0.4352307 0.3300142 0.4433984 0.4514112 0.4945122
#> [43] 0.4904250 0.5591064 0.5372852 0.3362153 0.5415618 0.5341986 0.3446184
#> [50] 0.5204829
# }
```
