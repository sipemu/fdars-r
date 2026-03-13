# Conformal Prediction for Logistic Regression

Conformal Prediction for Logistic Regression

## Usage

``` r
conformal.logistic(
  fdataobj,
  y,
  newdata,
  scalar.train = NULL,
  scalar.test = NULL,
  ncomp = 3,
  max.iter = 100,
  tol = 1e-06,
  score.type = c("lac", "aps"),
  cal.fraction = 0.25,
  alpha = 0.1,
  seed = NULL
)
```

## Arguments

- fdataobj:

  An object of class 'fdata' (training data).

- y:

  Binary response (0/1).

- newdata:

  An object of class 'fdata' (test data).

- scalar.train:

  Optional scalar covariates for training.

- scalar.test:

  Optional scalar covariates for test.

- ncomp:

  Number of FPC components (default 3).

- max.iter:

  Maximum iterations (default 100).

- tol:

  Convergence tolerance (default 1e-6).

- score.type:

  Nonconformity score: "lac" or "aps".

- cal.fraction:

  Calibration fraction (default 0.25).

- alpha:

  Miscoverage level (default 0.1).

- seed:

  Random seed.

## Value

Same as
[`conformal.classif`](https://sipemu.github.io/fdars-r/reference/conformal.classif.md).

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- rbinom(50, 1, 0.5)
cp <- conformal.logistic(fd[1:40, ], y[1:40], fd[41:50, ])
# }
```
