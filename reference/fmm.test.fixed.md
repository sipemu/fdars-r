# Permutation Test for Fixed Effects in FMM

Tests significance of fixed effects (covariates) using a
permutation-based F-test that respects the within-subject correlation.

## Usage

``` r
fmm.test.fixed(
  fdataobj,
  subject.ids,
  covariates,
  ncomp = 3,
  n.perm = 1000,
  seed = NULL
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- subject.ids:

  Integer vector of subject identifiers.

- covariates:

  Matrix of covariates (n_total x p).

- ncomp:

  Number of FPC components (default 3).

- n.perm:

  Number of permutations (default 1000).

- seed:

  Random seed.

## Value

An object of class 'fmm.test' with components:

- f.statistics:

  F-statistic per covariate

- p.values:

  P-values per covariate

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
subject <- rep(1:10, each = 5)
x <- cbind(rnorm(50))
test <- fmm.test.fixed(fd, subject.ids = subject, covariates = x,
                        n.perm = 100)
test
#> Fixed Effects Permutation Test
#> ==============================
#>   Permutations: 100 
#> 
#>             F.statistic P.value
#> Covariate 1      0.0029 0.80198
# }
```
