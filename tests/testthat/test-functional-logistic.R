# Tests for functional.logistic and predict.fregre.logistic

# =============================================================================
# Helper
# =============================================================================

.test_logistic_data <- function(n = 60, m = 20, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  for (i in 1:n) {
    X[i, ] <- sin(2 * pi * t) * (2 * (i > n / 2) - 1) + rnorm(m, sd = 0.3)
  }
  y <- as.numeric(seq_len(n) > n / 2)
  fd <- fdata(X, argvals = t)
  list(fd = fd, y = y)
}

# =============================================================================
# functional.logistic
# =============================================================================

test_that("functional.logistic returns correct structure", {
  d <- .test_logistic_data()
  fit <- functional.logistic(d$fd, d$y, ncomp = 3)

  expect_s3_class(fit, "fregre.logistic")
  expect_true(is.numeric(fit$intercept))
  expect_s3_class(fit$beta.t, "fdata")
  expect_equal(length(fit$probabilities), 60)
  expect_equal(length(fit$predicted.classes), 60)
  expect_true(fit$accuracy >= 0 && fit$accuracy <= 1)
  expect_true(is.finite(fit$log.likelihood))
  expect_true(is.finite(fit$aic))
  expect_true(is.finite(fit$bic))
  expect_equal(fit$ncomp, 3)
})

test_that("functional.logistic probabilities are in [0, 1]", {
  d <- .test_logistic_data()
  fit <- functional.logistic(d$fd, d$y, ncomp = 3)

  expect_true(all(fit$probabilities >= 0))
  expect_true(all(fit$probabilities <= 1))
})

test_that("functional.logistic predicted classes are 0 or 1", {
  d <- .test_logistic_data()
  fit <- functional.logistic(d$fd, d$y, ncomp = 3)

  expect_true(all(fit$predicted.classes %in% c(0, 1)))
})

test_that("functional.logistic achieves good accuracy on separable data", {
  d <- .test_logistic_data()
  fit <- functional.logistic(d$fd, d$y, ncomp = 3)

  expect_gt(fit$accuracy, 0.7)
})

test_that("functional.logistic with scalar covariates", {
  d <- .test_logistic_data()
  z <- matrix(rnorm(60 * 2), 60, 2)
  fit <- functional.logistic(d$fd, d$y, scalar.covariates = z, ncomp = 3)

  expect_s3_class(fit, "fregre.logistic")
  expect_true(length(fit$gamma) > 0)
})

test_that("functional.logistic validates input", {
  expect_error(functional.logistic("not_fdata", 1:10), "fdata")

  d <- .test_logistic_data()
  expect_error(functional.logistic(d$fd, 1:5), "Length of y")
})

test_that("functional.logistic stores FPCA internals", {
  d <- .test_logistic_data()
  fit <- functional.logistic(d$fd, d$y, ncomp = 3)

  expect_true(!is.null(fit$.fpca_mean))
  expect_true(!is.null(fit$.fpca_rotation))
  expect_true(!is.null(fit$.fpca_scores))
  expect_equal(ncol(fit$.fpca_scores), 3)
  expect_equal(nrow(fit$.fpca_scores), 60)
})

# =============================================================================
# predict.fregre.logistic
# =============================================================================

test_that("predict.fregre.logistic returns fitted probs when newdata=NULL", {
  d <- .test_logistic_data()
  fit <- functional.logistic(d$fd, d$y, ncomp = 3)

  pred <- predict(fit)
  expect_equal(pred, fit$probabilities)
})

test_that("predict.fregre.logistic with new data returns probs", {
  d <- .test_logistic_data()
  fit <- functional.logistic(d$fd, d$y, ncomp = 3)

  pred <- predict(fit, d$fd[1:10, ])
  expect_equal(length(pred), 10)
  expect_true(all(pred >= 0 & pred <= 1))
})

test_that("predict.fregre.logistic type='class' returns 0/1", {
  d <- .test_logistic_data()
  fit <- functional.logistic(d$fd, d$y, ncomp = 3)

  classes <- predict(fit, d$fd[1:10, ], type = "class")
  expect_equal(length(classes), 10)
  expect_true(all(classes %in% c(0, 1)))
})

test_that("predict.fregre.logistic with scalar covariates", {
  d <- .test_logistic_data()
  z <- matrix(rnorm(60 * 2), 60, 2)
  fit <- functional.logistic(d$fd, d$y, scalar.covariates = z, ncomp = 3)

  pred <- predict(fit, d$fd[1:10, ], new.scalar = z[1:10, ])
  expect_equal(length(pred), 10)
  expect_true(all(pred >= 0 & pred <= 1))
})

test_that("predict.fregre.logistic converts matrix to fdata", {
  d <- .test_logistic_data()
  fit <- functional.logistic(d$fd, d$y, ncomp = 3)

  # Pass raw matrix instead of fdata
  pred <- predict(fit, d$fd$data[1:5, , drop = FALSE])
  expect_equal(length(pred), 5)
  expect_true(all(pred >= 0 & pred <= 1))
})

# =============================================================================
# print.fregre.logistic
# =============================================================================

test_that("print.fregre.logistic works", {
  d <- .test_logistic_data()
  fit <- functional.logistic(d$fd, d$y, ncomp = 3)

  expect_output(print(fit), "Logistic")
})
