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
#>  [1]  9.157067  8.691053 10.459906 12.457908 14.401164 14.187978 15.723110
#>  [8] 17.617270 18.973128 20.584498 20.584498 20.584498 20.584498 20.584498
#> [15] 20.584498
#> 
#> $bic
#>  [1] 12.98111 14.42712 18.10800 22.01802 25.87330 27.57214 31.01929 34.82548
#>  [9] 38.09336 41.61675 41.61675 41.61675 41.61675 41.61675 41.61675
#> 
#> $gcv
#>  [1] 1.199439 1.183961 1.207095 1.254431 1.330111 1.324644 1.392962 1.507716
#>  [9] 1.556828 1.635819 1.635819 1.635819 1.635819 1.635819 1.635819
#> 
# }
```
