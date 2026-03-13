# Cross-Conformal (CV+) Classification

Constructs conformal prediction sets using K-fold cross-conformal
inference.

## Usage

``` r
cv.conformal.classification(
  fdataobj,
  y,
  newdata,
  covariates.train = NULL,
  covariates.test = NULL,
  ncomp = 3,
  classifier = "lda",
  k.nn = 5,
  score.type = c("lac", "aps"),
  n.folds = 5,
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

- n.folds:

  Number of folds (default 5).

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
cp <- cv.conformal.classification(fd[1:40, ], y[1:40], fd[41:50, ])
# }
```
