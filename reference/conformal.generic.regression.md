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

- cal.fraction:

  Fraction of data for calibration (default 0.25).

- alpha:

  Miscoverage level (default 0.1).

- seed:

  Random seed.

## Value

Same as
[`conformal.fregre.lm`](https://sipemu.github.io/fdars-r/reference/conformal.fregre.lm.md).

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
model <- fregre.lm(fd, y)
cp <- conformal.generic.regression(model, fd, y, fd[1:10, ])
# }
```
