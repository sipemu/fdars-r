# Functional ANOVA

Tests for differences in mean functions across groups using a
permutation-based F-test.

## Usage

``` r
fanova(fdataobj, groups, n.perm = 1000)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- groups:

  Integer vector of group labels.

- n.perm:

  Number of permutations (default 1000).

## Value

An object of class 'fanova' with components:

- group.means:

  Group mean functions as fdata

- overall.mean:

  Overall mean function

- f.statistic.t:

  Pointwise F-statistic

- global.statistic:

  Global test statistic

- p.value:

  P-value from permutation test

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
groups <- rep(1:2, each = 25)
res <- fanova(fd, groups = groups, n.perm = 100)
res
#> Functional ANOVA
#> ================
#>   Number of groups: 2 
#>   Number of observations: 50 
#>   Global F-statistic: 0.8663 
#>   P-value: 0.75248 
#>   Permutations: 100 
# }
```
