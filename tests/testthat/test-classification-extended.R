# Extended tests for classification: SVM, covariates, edge cases

# =============================================================================
# Helper
# =============================================================================

.test_classif_data <- function(n_per = 20, m = 20, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X1 <- matrix(0, n_per, m)
  X2 <- matrix(0, n_per, m)
  for (i in 1:n_per) {
    X1[i, ] <- sin(2 * pi * t) + rnorm(m, sd = 0.2)
    X2[i, ] <- cos(2 * pi * t) + rnorm(m, sd = 0.2)
  }
  X <- rbind(X1, X2)
  y <- rep(1:2, each = n_per)
  fd <- fdata(X, argvals = t)
  list(fd = fd, y = y, n = 2 * n_per)
}

# =============================================================================
# fclassif with SVM
# =============================================================================

test_that("fclassif with SVM returns correct structure", {
  skip_if_not_installed("e1071")
  d <- .test_classif_data()
  result <- fclassif(d$fd, d$y, method = "svm", ncomp = 3)

  expect_s3_class(result, "fclassif")
  expect_equal(result$method, "svm")
  expect_equal(length(result$predicted), d$n)
  expect_true(all(result$predicted %in% 1:2))
  expect_true(result$accuracy >= 0 && result$accuracy <= 1)
})

test_that("fclassif SVM achieves reasonable accuracy on separable data", {
  skip_if_not_installed("e1071")
  d <- .test_classif_data()
  result <- fclassif(d$fd, d$y, method = "svm", ncomp = 3)

  expect_gt(result$accuracy, 0.7)
})

test_that("fclassif SVM confusion matrix is valid", {
  skip_if_not_installed("e1071")
  d <- .test_classif_data()
  result <- fclassif(d$fd, d$y, method = "svm", ncomp = 3)

  expect_true(is.matrix(result$confusion))
  expect_equal(sum(result$confusion), d$n)
})

test_that("fclassif.cv with SVM returns correct structure", {
  skip_if_not_installed("e1071")
  d <- .test_classif_data()
  result <- fclassif.cv(d$fd, d$y, method = "svm", ncomp = 3,
                         nfold = 3, seed = 42)

  expect_s3_class(result, "fclassif.cv")
  expect_true("error.rate" %in% names(result))
  expect_gte(result$error.rate, 0)
  expect_lte(result$error.rate, 1)
})

# =============================================================================
# fclassif with 3 classes
# =============================================================================

test_that("fclassif handles 3-class problems", {
  set.seed(42)
  m <- 20
  t <- seq(0, 1, length.out = m)
  n_per <- 15
  X <- matrix(0, 3 * n_per, m)
  for (i in 1:n_per) {
    X[i, ] <- sin(2 * pi * t) + rnorm(m, sd = 0.15)
    X[n_per + i, ] <- cos(2 * pi * t) + rnorm(m, sd = 0.15)
    X[2 * n_per + i, ] <- t + rnorm(m, sd = 0.15)
  }
  fd <- fdata(X, argvals = t)
  y <- rep(1:3, each = n_per)

  result_lda <- fclassif(fd, y, method = "lda", ncomp = 3)
  expect_equal(length(unique(result_lda$predicted)), 3)
  expect_equal(nrow(result_lda$confusion), 3)
})

# =============================================================================
# fclassif with covariates
# =============================================================================

test_that("fclassif with covariates works for LDA", {
  d <- .test_classif_data()
  z <- matrix(rnorm(d$n * 2), d$n, 2)
  result <- fclassif(d$fd, d$y, method = "lda", covariates = z, ncomp = 3)

  expect_s3_class(result, "fclassif")
  expect_equal(length(result$predicted), d$n)
})

test_that("fclassif with covariates works for kNN", {
  d <- .test_classif_data()
  z <- matrix(rnorm(d$n * 2), d$n, 2)
  result <- fclassif(d$fd, d$y, method = "knn", covariates = z, ncomp = 3, k = 5)

  expect_s3_class(result, "fclassif")
  expect_equal(length(result$predicted), d$n)
})

# =============================================================================
# fclassif.cv edge cases
# =============================================================================

test_that("fclassif.cv with different nfold values", {
  d <- .test_classif_data()

  cv3 <- fclassif.cv(d$fd, d$y, method = "lda", ncomp = 3, nfold = 3, seed = 42)
  cv5 <- fclassif.cv(d$fd, d$y, method = "lda", ncomp = 3, nfold = 5, seed = 42)

  expect_s3_class(cv3, "fclassif.cv")
  expect_s3_class(cv5, "fclassif.cv")
  expect_true(length(cv3$fold.errors) == 3 || is.numeric(cv3$error.rate))
  expect_true(length(cv5$fold.errors) == 5 || is.numeric(cv5$error.rate))
})

test_that("fclassif.cv with kNN method", {
  d <- .test_classif_data()
  cv <- fclassif.cv(d$fd, d$y, method = "knn", ncomp = 3, nfold = 3,
                     seed = 42, k = 5)

  expect_s3_class(cv, "fclassif.cv")
  expect_gte(cv$error.rate, 0)
  expect_lte(cv$error.rate, 1)
})

test_that("fclassif.cv with QDA method", {
  d <- .test_classif_data()
  cv <- fclassif.cv(d$fd, d$y, method = "qda", ncomp = 3, nfold = 3, seed = 42)

  expect_s3_class(cv, "fclassif.cv")
  expect_gte(cv$error.rate, 0)
  expect_lte(cv$error.rate, 1)
})
