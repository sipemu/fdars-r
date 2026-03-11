# Elastic Principal Component Regression

Fits scalar-on-function regression using elastic FPCA scores as
predictors.

## Usage

``` r
elastic.pcr(
  fdataobj,
  y,
  ncomp = 3,
  pca.method = c("vertical", "horizontal", "joint"),
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

- ncomp:

  Number of principal components (default 3).

- pca.method:

  PCA decomposition method: "vertical" (amplitude), "horizontal"
  (phase), or "joint" (combined).

- lambda:

  Regularization parameter for alignment (default 0).

- max.iter:

  Maximum alignment iterations (default 20).

- tol:

  Convergence tolerance (default 1e-4).

## Value

An object of class 'elastic.pcr' with components:

- alpha:

  Intercept

- coefficients:

  PC regression coefficients

- fitted.values:

  Fitted response values

- sse:

  Sum of squared errors

- r.squared:

  R-squared

- pca.method:

  PCA method used

- karcher.mean:

  Karcher mean curve

- vert.scores:

  Vertical FPCA scores (if applicable)

- horiz.scores:

  Horizontal FPCA scores (if applicable)
