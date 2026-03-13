# Generic Conformal Classification

Split-conformal prediction sets using a pre-fitted functional.logistic
model.

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

  Binary response (0/1).

- newdata:

  An object of class 'fdata' (test data).

- scalar.train:

  Optional scalar covariates for training.

- scalar.test:

  Optional scalar covariates for test.

- score.type:

  Nonconformity score: "lac" or "aps" (default "lac").

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
model <- functional.logistic(fd, y)
cp <- conformal.generic.classification(model, fd, y, fd[1:10, ])
# }
```
