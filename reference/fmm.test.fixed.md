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
