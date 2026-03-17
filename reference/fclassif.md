# Supervised Classification of Functional Data

Functions for classifying functional observations into discrete groups.
Functional Classification

## Usage

``` r
fclassif(
  fdataobj,
  y,
  method = c("lda", "qda", "knn", "kernel", "dd", "svm"),
  covariates = NULL,
  ncomp = 3,
  ...
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- y:

  Integer vector of class labels (1-indexed).

- method:

  Classification method: "lda", "qda", "knn", "kernel", "dd", or "svm".

- covariates:

  Optional matrix of scalar covariates.

- ncomp:

  Number of FPC components (default 3). Used by lda, qda, knn, svm.

- ...:

  Additional arguments:

  k

  :   Number of neighbors for kNN (default 5).

  h.func

  :   Bandwidth for functional kernel (default auto).

  h.scalar

  :   Bandwidth for scalar kernel (default auto).

  kernel

  :   SVM kernel type: "radial" (default), "linear", "polynomial",
      "sigmoid".

  cost

  :   SVM cost parameter (default 1).

  gamma

  :   SVM kernel parameter (default 1/ncomp).

## Value

An object of class 'fclassif' with components:

- predicted:

  Predicted class labels

- probabilities:

  Posterior probabilities (if available)

- accuracy:

  Training accuracy

- confusion:

  Confusion matrix

- method:

  Method used

- ncomp:

  Number of FPC components

## Details

Classifies functional data using one of several methods: LDA, QDA, kNN,
kernel, DD-plot (depth-vs-depth), or SVM (support vector machine).

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
y <- rep(1:2, each = 25)
fit <- fclassif(fd, y, method = "lda", ncomp = 3)
fit$accuracy
#> [1] 0.58
# }
```
