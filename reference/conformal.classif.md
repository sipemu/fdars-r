# Conformal Classification

Constructs conformal prediction sets for functional classification.

## Usage

``` r
conformal.classif(
  fdataobj,
  y,
  newdata,
  covariates.train = NULL,
  covariates.test = NULL,
  ncomp = 3,
  classifier = "lda",
  k.nn = 5,
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

  Class labels (integer vector, 0-indexed).

- newdata:

  An object of class 'fdata' (test data).

- covariates.train:

  Optional scalar covariates for training.

- covariates.test:

  Optional scalar covariates for test.

- ncomp:

  Number of FPC components (default 3).

- classifier:

  Classifier: "lda", "qda", or "knn" (default "lda").

- k.nn:

  k for kNN classifier (default 5).

- score.type:

  Nonconformity score: "lac" or "aps" (default "lac").

- cal.fraction:

  Calibration fraction (default 0.25).

- alpha:

  Miscoverage level (default 0.1).

- seed:

  Random seed.

## Value

A list with components:

- predicted.classes:

  Point predictions

- set.sizes:

  Size of each prediction set

- average.set.size:

  Mean prediction set size

- coverage:

  Empirical coverage

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- sample(0:1, 50, replace = TRUE)
cp <- conformal.classif(fd[1:40, ], y[1:40], fd[41:50, ])
# }
```
