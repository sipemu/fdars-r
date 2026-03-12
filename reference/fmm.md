# Functional Mixed Models

Functions for fitting functional mixed models with subject-level random
effects. Functional Mixed Model

## Usage

``` r
fmm(fdataobj, subject.ids, covariates = NULL, ncomp = 3)
```

## Arguments

- fdataobj:

  An object of class 'fdata'. All observations stacked (multiple per
  subject).

- subject.ids:

  Integer vector of subject identifiers (same length as nrow(data)).

- covariates:

  Optional matrix of covariates (n_total x p).

- ncomp:

  Number of FPC components for random effects (default 3).

## Value

An object of class 'fmm' with components:

- mean.function:

  Overall mean function as fdata

- beta.functions:

  Fixed effect coefficient functions as fdata

- random.effects:

  Random effect functions per subject as fdata

- fitted:

  Fitted values as fdata

- residuals:

  Residuals as fdata

- sigma2.eps:

  Residual variance estimate

- sigma2.u:

  Random effect variance estimates

- n.subjects:

  Number of unique subjects

## Details

Fits a functional mixed model for repeated measures data:
`Y_ij(t) = mu(t) + sum_k x_ijk beta_k(t) + b_i(t) + eps_ij(t)` where
b_i(t) are subject-level random effects.

## Examples

``` r
# \donttest{
# 10 subjects, 5 curves each = 50 total curves
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
subject <- rep(1:10, each = 5)
fit <- fmm(fd, subject.ids = subject)
fit
#> Functional Mixed Model
#> ======================
#>   Number of observations: 50 
#>   Number of subjects: 10 
#>   FPC components: 3 
#>   Residual variance (sigma2_eps): 0.170716 
# }
```
