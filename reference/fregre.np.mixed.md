# Nonparametric Functional Regression with Mixed Predictors

Fits a nonparametric kernel regression with functional and scalar
predictors using the Rust backend.

## Usage

``` r
fregre.np.mixed(
  fdataobj,
  y,
  scalar.covariates = NULL,
  h.func = NULL,
  h.scalar = NULL
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- y:

  Response vector.

- scalar.covariates:

  Optional matrix of scalar covariates.

- h.func:

  Bandwidth for functional distance kernel. If NULL, auto-selected.

- h.scalar:

  Bandwidth for scalar covariates kernel. If NULL, auto-selected.

## Value

A fitted regression object of class 'fregre.np'.
