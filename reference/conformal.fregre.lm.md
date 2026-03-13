# Conformal Prediction for Functional Linear Model

Constructs prediction intervals using split conformal inference for the
FPC-based functional linear model.

## Usage

``` r
conformal.fregre.lm(
  fdataobj,
  y,
  newdata,
  scalar.train = NULL,
  scalar.test = NULL,
  ncomp = 3,
  cal.fraction = 0.25,
  alpha = 0.1,
  seed = NULL
)
```

## Arguments

- fdataobj:

  An object of class 'fdata' (training data).

- y:

  Response vector (training).

- newdata:

  An object of class 'fdata' (test data).

- scalar.train:

  Optional scalar covariates for training.

- scalar.test:

  Optional scalar covariates for test.

- ncomp:

  Number of FPC components (default 3).

- cal.fraction:

  Fraction of data for calibration (default 0.25).

- alpha:

  Miscoverage level (default 0.1 for 90 percent intervals).

- seed:

  Random seed.

## Value

A list with components:

- predictions:

  Point predictions for test data

- lower:

  Lower bounds of prediction intervals

- upper:

  Upper bounds of prediction intervals

- residual.quantile:

  Calibration residual quantile

- coverage:

  Empirical coverage on calibration set

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
cp <- conformal.fregre.lm(fd[1:40, ], y[1:40], fd[41:50, ], ncomp = 3)
# }
```
