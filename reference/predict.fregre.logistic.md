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
#>  [1] 0.2315648 0.4543441 0.2208339 0.4758234 0.3836514 0.4386480 0.3852055
#>  [8] 0.4941099 0.4766884 0.3229255 0.4572099 0.3375068 0.3679262 0.6786208
#> [15] 0.3323573 0.4638493 0.3731966 0.4990061 0.5790934 0.4482920 0.3040211
#> [22] 0.6046318 0.4643096 0.5931189 0.5611110 0.4581105 0.4738770 0.3905563
#> [29] 0.2501936 0.5781530 0.4509228 0.3678800 0.4109633 0.3909392 0.5444729
#> [36] 0.3860902 0.2425125 0.4245526 0.3538779 0.4236073 0.4257720 0.5066651
#> [43] 0.4997009 0.6493563 0.4917168 0.3424366 0.5089224 0.5193545 0.4241595
#> [50] 0.5371610
# }
```
