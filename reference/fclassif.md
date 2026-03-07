# Supervised Classification of Functional Data

Functions for classifying functional observations into discrete groups.
Functional Classification

## Usage

``` r
fclassif(
  fdataobj,
  y,
  method = c("lda", "qda", "knn", "kernel", "dd"),
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

  Classification method: "lda", "qda", "knn", "kernel", or "dd".

- covariates:

  Optional matrix of scalar covariates.

- ncomp:

  Number of FPC components (default 3). Used by lda, qda, knn.

- ...:

  Additional arguments:

  k

  :   Number of neighbors for kNN (default 5).

  h.func

  :   Bandwidth for functional kernel (default auto).

  h.scalar

  :   Bandwidth for scalar kernel (default auto).

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
kernel, or DD-plot (depth-vs-depth).
