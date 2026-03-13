# Conformal Prediction for Nonparametric Functional Regression

Conformal Prediction for Nonparametric Functional Regression

## Usage

``` r
conformal.fregre.np(
  fdataobj,
  y,
  newdata,
  scalar.train = NULL,
  scalar.test = NULL,
  h.func = 0,
  h.scalar = 0,
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

- h.func:

  Functional bandwidth (0 for auto).

- h.scalar:

  Scalar bandwidth (0 for auto).

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
cp <- conformal.fregre.np(fd[1:40, ], y[1:40], fd[41:50, ])
# }
```
