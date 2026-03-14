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
#>  [1]  9.234859 11.053008 12.944234 13.774887 14.589334 16.062683 16.952139
#>  [8] 18.164922 20.134115 20.584498 20.584498 20.584498 20.584498 20.584498
#> [15] 20.584498
#> 
#> $bic
#>  [1] 13.05890 16.78908 20.59233 23.33500 26.06147 29.44684 32.24832 35.37313
#>  [9] 39.25435 41.61675 41.61675 41.61675 41.61675 41.61675 41.61675
#> 
#> $gcv
#>  [1] 1.183827 1.250213 1.278733 1.281603 1.326597 1.420450 1.444756 1.491719
#>  [9] 1.569823 1.635819 1.635819 1.635819 1.635819 1.635819 1.635819
#> 
# }
```
