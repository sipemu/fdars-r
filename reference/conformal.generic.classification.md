# Generic Conformal Classification

Split-conformal prediction sets using a pre-fitted functional.logistic
model. Only supports binary classification (2 classes).

## Usage

``` r
conformal.generic.classification(
  model,
  fdataobj,
  y,
  newdata,
  scalar.train = NULL,
  scalar.test = NULL,
  score.type = c("lac", "aps"),
  calibration.indices = NULL,
  cal.fraction = 0.25,
  alpha = 0.1,
  seed = NULL
)
```

## Arguments

- model:

  A fitted model object from
  [`functional.logistic`](https://sipemu.github.io/fdars-r/reference/functional.logistic.md).

- fdataobj:

  An object of class 'fdata' (training data).

- y:

  Binary response (0/1). Only binary (2-class) classification is
  supported; for multiclass use
  [`conformal.classif`](https://sipemu.github.io/fdars-r/reference/conformal.classif.md).

- newdata:

  An object of class 'fdata' (test data).

- scalar.train:

  Optional scalar covariates for training.

- scalar.test:

  Optional scalar covariates for test.

- score.type:

  Nonconformity score: "lac" or "aps" (default "lac").

- calibration.indices:

  Optional integer vector of 1-based indices into the training data to
  use as the calibration set. When provided, these observations should
  have been held out during model fitting so that calibration scores are
  out-of-sample, restoring the coverage guarantee.

- cal.fraction:

  Calibration fraction (default 0.25). Ignored when
  `calibration.indices` is provided.

- alpha:

  Miscoverage level (default 0.1).

- seed:

  Random seed.

## Value

Same as
[`conformal.classif`](https://sipemu.github.io/fdars-r/reference/conformal.classif.md).

## Warning

The model was trained on ALL data including the calibration subset, so
calibration scores are in-sample and the coverage guarantee is broken.
Use
[`conformal.logistic`](https://sipemu.github.io/fdars-r/reference/conformal.logistic.md)
or
[`cv.conformal.classification`](https://sipemu.github.io/fdars-r/reference/cv.conformal.classification.md)
for valid coverage. For multiclass problems, use
[`conformal.classif`](https://sipemu.github.io/fdars-r/reference/conformal.classif.md).

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- rbinom(50, 1, 0.5)
model <- functional.logistic(fd, y)
cp <- conformal.generic.classification(model, fd, y, fd[1:10, ])
#> Warning: conformal.generic.classification uses the pre-fitted model without refitting. Calibration scores are in-sample, so coverage guarantee is broken. Supply calibration.indices (held-out indices) for valid coverage, or use conformal.logistic() / cv.conformal.classification() instead.
# }
```
