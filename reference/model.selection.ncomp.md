# Model Selection for Number of FPC Components

Selects the optimal number of FPC components for scalar-on-function
regression using AIC, BIC, or GCV criterion.

## Usage

``` r
model.selection.ncomp(
  fdataobj,
  y,
  scalar.covariates = NULL,
  max.ncomp = 15,
  criterion = c("aic", "bic", "gcv")
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- y:

  Response vector (scalar).

- scalar.covariates:

  Optional matrix of scalar covariates.

- max.ncomp:

  Maximum number of components to evaluate (default 15).

- criterion:

  Selection criterion: "aic", "bic", or "gcv".

## Value

A list with components:

- best.ncomp:

  Optimal number of components

- ncomp:

  Vector of tested component counts

- aic:

  AIC values

- bic:

  BIC values

- gcv:

  GCV values

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
model.selection.ncomp(fd, y, criterion = "bic")
#> $best.ncomp
#> [1] 1
#> 
#> $ncomp
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
#> 
#> $aic
#>  [1]  3.164190  5.104452  7.104450  9.073390 11.027029 12.714502 14.711688
#>  [8] 12.583049 13.117048 15.114595 15.114595 15.114595 15.114595 15.114595
#> [15] 15.114595
#> 
#> $bic
#>  [1]  6.988236 10.840521 14.752542 18.633505 22.499167 26.098663 30.007872
#>  [8] 29.791256 32.237278 36.146848 36.146848 36.146848 36.146848 36.146848
#> [15] 36.146848
#> 
#> $gcv
#>  [1] 1.079935 1.119478 1.182352 1.210059 1.265301 1.309340 1.368977 1.322212
#>  [9] 1.324895 1.382475 1.382475 1.382475 1.382475 1.382475 1.382475
#> 
# }
```
