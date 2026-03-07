# Tests for functional classification

# =============================================================================
# Test Data Setup
# =============================================================================

create_classif_data <- function(n_per_class = 15, m = 30, seed = 42) {
  set.seed(seed)
  t_grid <- seq(0, 1, length.out = m)

  # Class 1: sine-like curves
  data1 <- matrix(0, nrow = n_per_class, ncol = m)
  for (i in 1:n_per_class) {
    data1[i, ] <- sin(2 * pi * t_grid) + rnorm(m, 0, 0.15)
  }

  # Class 2: cosine-like curves
  data2 <- matrix(0, nrow = n_per_class, ncol = m)
  for (i in 1:n_per_class) {
    data2[i, ] <- cos(2 * pi * t_grid) + rnorm(m, 0, 0.15)
  }

  # Class 3: linear curves
  data3 <- matrix(0, nrow = n_per_class, ncol = m)
  for (i in 1:n_per_class) {
    data3[i, ] <- t_grid + rnorm(m, 0, 0.15)
  }

  data_mat <- rbind(data1, data2, data3)
  y <- rep(1:3, each = n_per_class)

  list(
    fd = fdars::fdata(data_mat, argvals = t_grid),
    y = y,
    n = 3 * n_per_class,
    m = m
  )
}

# =============================================================================
# fclassif Tests
# =============================================================================

test_that("fclassif with LDA returns correct structure", {
  data <- create_classif_data()
  result <- fdars::fclassif(data$fd, data$y, method = "lda", ncomp = 3)

  expect_s3_class(result, "fclassif")
  expect_true("predicted" %in% names(result))
  expect_true("accuracy" %in% names(result))
  expect_true("confusion" %in% names(result))
  expect_true("method" %in% names(result))
  expect_equal(result$method, "lda")
  expect_length(result$predicted, data$n)
})

test_that("fclassif with QDA returns correct structure", {
  data <- create_classif_data()
  result <- fdars::fclassif(data$fd, data$y, method = "qda", ncomp = 3)

  expect_s3_class(result, "fclassif")
  expect_equal(result$method, "qda")
  expect_length(result$predicted, data$n)
})

test_that("fclassif with kNN returns correct structure", {
  data <- create_classif_data()
  result <- fdars::fclassif(data$fd, data$y, method = "knn", ncomp = 3, k = 5)

  expect_s3_class(result, "fclassif")
  expect_equal(result$method, "knn")
  expect_length(result$predicted, data$n)
})

test_that("fclassif with kernel returns correct structure", {
  data <- create_classif_data()
  result <- fdars::fclassif(data$fd, data$y, method = "kernel")

  expect_s3_class(result, "fclassif")
  expect_equal(result$method, "kernel")
  expect_length(result$predicted, data$n)
})

test_that("fclassif with DD returns correct structure", {
  data <- create_classif_data()
  result <- fdars::fclassif(data$fd, data$y, method = "dd")

  expect_s3_class(result, "fclassif")
  expect_equal(result$method, "dd")
  expect_length(result$predicted, data$n)
})

test_that("fclassif achieves reasonable accuracy on separable data", {
  data <- create_classif_data()
  result <- fdars::fclassif(data$fd, data$y, method = "lda", ncomp = 3)

  expect_gt(result$accuracy, 0.5)
})

test_that("fclassif confusion matrix is valid", {
  data <- create_classif_data()
  result <- fdars::fclassif(data$fd, data$y, method = "lda", ncomp = 3)

  conf <- result$confusion
  expect_true(is.matrix(conf))
  expect_equal(sum(conf), data$n)
})

# =============================================================================
# fclassif.cv Tests
# =============================================================================

test_that("fclassif.cv returns correct structure", {
  data <- create_classif_data()
  result <- fdars::fclassif.cv(data$fd, data$y, method = "lda",
                                ncomp = 3, nfold = 5, seed = 123)

  expect_s3_class(result, "fclassif.cv")
  expect_true("error.rate" %in% names(result))
  expect_true("fold.errors" %in% names(result))
  expect_true("best.ncomp" %in% names(result))
})

test_that("fclassif.cv error rate is valid", {
  data <- create_classif_data()
  result <- fdars::fclassif.cv(data$fd, data$y, method = "lda",
                                ncomp = 3, nfold = 5, seed = 123)

  expect_gte(result$error.rate, 0)
  expect_lte(result$error.rate, 1)
})

# =============================================================================
# Input Validation Tests
# =============================================================================

test_that("fclassif rejects non-fdata input", {
  X <- matrix(rnorm(100), 10, 10)
  expect_error(fdars::fclassif(X, rep(1:2, 5)))
})

test_that("fclassif rejects mismatched y length", {
  fd <- fdars::fdata(matrix(rnorm(100), 10, 10))
  expect_error(fdars::fclassif(fd, 1:5, method = "lda"), "Length of y")
})

# =============================================================================
# Print and Plot Tests
# =============================================================================

test_that("print.fclassif produces output", {
  data <- create_classif_data()
  result <- fdars::fclassif(data$fd, data$y, method = "lda", ncomp = 3)

  expect_output(print(result), "Classification")
  expect_output(print(result), "LDA")
})

test_that("print.fclassif.cv produces output", {
  data <- create_classif_data()
  result <- fdars::fclassif.cv(data$fd, data$y, method = "lda",
                                ncomp = 3, nfold = 5, seed = 123)

  expect_output(print(result), "Cross-Validated")
  expect_output(print(result), "Error rate")
})

test_that("plot.fclassif works", {
  data <- create_classif_data()
  result <- fdars::fclassif(data$fd, data$y, method = "lda", ncomp = 3)

  expect_no_error(plot(result))
})
