# Generic Conformal Regression

Split-conformal prediction intervals using a pre-fitted fregre.lm model.
Uses the model's FPCA components for fast prediction without refitting.

## Usage

``` r
conformal.generic.regression(
  model,
  fdataobj,
  y,
  newdata,
  scalar.train = NULL,
  scalar.test = NULL,
  calibration.indices = NULL,
  cal.fraction = 0.25,
  alpha = 0.1,
  seed = NULL
)
```

## Arguments

- model:

  A fitted `fregre.lm` model object.

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

- calibration.indices:

  Optional integer vector of 1-based indices into the training data to
  use as the calibration set. When provided, these observations should
  have been held out during model fitting so that calibration residuals
  are out-of-sample, restoring the coverage guarantee. If NULL
  (default), calibration indices are randomly selected from all training
  data (in-sample).

- cal.fraction:

  Fraction of data for calibration (default 0.25). Ignored when
  `calibration.indices` is provided.

- alpha:

  Miscoverage level (default 0.1).

- seed:

  Random seed.

## Value

Same as
[`conformal.fregre.lm`](https://sipemu.github.io/fdars-r/reference/conformal.fregre.lm.md).

## Warning

The model was trained on ALL data including the calibration subset, so
calibration residuals are in-sample and systematically too small. The
distribution-free coverage guarantee is broken. Use
[`conformal.fregre.lm`](https://sipemu.github.io/fdars-r/reference/conformal.fregre.lm.md)
or
[`cv.conformal.regression`](https://sipemu.github.io/fdars-r/reference/cv.conformal.regression.md)
for valid coverage.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
model <- fregre.lm(fd, y)
cp <- conformal.generic.regression(model, fd, y, fd[1:10, ])
#> Warning: conformal.generic.regression uses the pre-fitted model without refitting. Calibration residuals are in-sample, so coverage guarantee is broken. Supply calibration.indices (held-out indices) for valid coverage, or use conformal.fregre.lm() / cv.conformal.regression() instead.
# }
```
