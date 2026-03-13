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
#>  [1] 0.3672389 0.1906722 0.4567561 0.3944575 0.4256876 0.4597095 0.4590716
#>  [8] 0.5900596 0.3340129 0.3232835 0.4854059 0.3451131 0.5181361 0.4607318
#> [15] 0.3888734 0.4434599 0.3741979 0.5287117 0.4159629 0.3170818 0.4091006
#> [22] 0.4321281 0.6921541 0.4434681 0.3033496 0.4175073 0.3907801 0.5265094
#> [29] 0.6050095 0.3418941 0.4099263 0.4482277 0.5103996 0.4072364 0.5382773
#> [36] 0.2280286 0.3112116 0.3454929 0.4332496 0.3678997 0.5337529 0.3592158
#> [43] 0.3841660 0.3150101 0.5096921 0.2915168 0.7049031 0.3515220 0.3900477
#> [50] 0.3196992
# }
```
