# Elastic Regression

Alignment-integrated regression and classification for functional data.
Elastic Scalar-on-Function Regression

## Usage

``` r
elastic.regression(
  fdataobj,
  y,
  ncomp.beta = 10,
  lambda = 0,
  max.iter = 20,
  tol = 1e-04
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- y:

  Response vector (numeric).

- ncomp.beta:

  Number of basis functions for the beta coefficient (default 10).

- lambda:

  Regularization parameter (default 0).

- max.iter:

  Maximum iterations for the alignment-regression loop (default 20).

- tol:

  Convergence tolerance (default 1e-4).

## Value

An object of class 'elastic.regression' with components:

- alpha:

  Intercept

- beta:

  Beta coefficient function

- fitted.values:

  Fitted response values

- residuals:

  Residuals

- sse:

  Sum of squared errors

- r.squared:

  R-squared

- gammas:

  Estimated warping functions

- aligned.srsfs:

  Aligned SRSF transforms

- n.iter:

  Number of iterations to convergence

- fdataobj:

  Original functional data

- y:

  Response vector

## Details

Fits a functional linear model with simultaneous elastic alignment of
the functional covariate. The alignment and regression are jointly
optimized.
