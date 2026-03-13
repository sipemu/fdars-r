# Conformal Prediction for Elastic Regression

Conformal Prediction for Elastic Regression

## Usage

``` r
conformal.elastic.regression(
  fdataobj,
  y,
  newdata,
  ncomp.beta = 10,
  lambda = 0,
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

- ncomp.beta:

  Number of basis functions for beta (default 10).

- lambda:

  Regularization (default 0).

- cal.fraction:

  Calibration fraction (default 0.25).

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
cp <- conformal.elastic.regression(fd[1:40, ], y[1:40], fd[41:50, ])
# }
```
