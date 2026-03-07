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

Fits a functional mixed model for repeated measures data: Y_ij(t) =
mu(t) + sum_k x_ij,k beta_k(t) + b_i(t) + eps_ij(t) where b_i(t) are
subject-level random effects.
