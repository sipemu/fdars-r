# Conformal Prediction for Elastic Logistic Regression

Conformal Prediction for Elastic Logistic Regression

## Usage

``` r
conformal.elastic.logistic(
  fdataobj,
  y,
  newdata,
  lambda = 0,
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

- lambda:

  Regularization parameter (default 0).

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
y <- sample(0:1, 50, replace = TRUE)
cp <- conformal.elastic.logistic(fd[1:40, ], y[1:40], fd[41:50, ])
# }
```
