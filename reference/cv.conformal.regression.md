# Cross-Conformal (CV+) Regression

Constructs prediction intervals using K-fold cross-conformal inference.
No data is wasted on a calibration split; each fold provides out-of-fold
calibration residuals.

## Usage

``` r
cv.conformal.regression(
  fdataobj,
  y,
  newdata,
  method = c("fregre.lm", "fregre.np"),
  scalar.train = NULL,
  scalar.test = NULL,
  ncomp = 3,
  h.func = 0,
  h.scalar = 0,
  n.folds = 5,
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

- method:

  Regression method: "fregre.lm" or "fregre.np".

- scalar.train:

  Optional scalar covariates for training.

- scalar.test:

  Optional scalar covariates for test.

- ncomp:

  Number of FPC components (default 3, for fregre.lm).

- h.func:

  Functional bandwidth (default 0 = auto, for fregre.np).

- h.scalar:

  Scalar bandwidth (default 0 = auto, for fregre.np).

- n.folds:

  Number of folds (default 5).

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
cp <- cv.conformal.regression(fd[1:40, ], y[1:40], fd[41:50, ])
# }
```
