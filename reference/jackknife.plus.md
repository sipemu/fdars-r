# Jackknife+ Regression

Constructs prediction intervals using jackknife+ (leave-one-out
conformal). Most sample-efficient but most expensive method (n refits).

## Usage

``` r
jackknife.plus(
  fdataobj,
  y,
  newdata,
  method = c("fregre.lm", "fregre.np"),
  scalar.train = NULL,
  scalar.test = NULL,
  ncomp = 3,
  h.func = 0,
  h.scalar = 0,
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

- alpha:

  Miscoverage level (default 0.1 for 90 percent intervals).

- seed:

  Random seed (unused, kept for API consistency).

## Value

Same as
[`cv.conformal.regression`](https://sipemu.github.io/fdars-r/reference/cv.conformal.regression.md).

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
jp <- jackknife.plus(fd[1:40, ], y[1:40], fd[41:50, ])
# }
```
