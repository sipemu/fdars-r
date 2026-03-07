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
