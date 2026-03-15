# Tests for robust functional regression (fregre.l1, fregre.huber)

# Helper: create test data with a clear linear signal
.test_robust_data <- function(n = 40, m = 20, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1, 1, 0.3)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(m, sd = 0.1)
    y[i] <- 2 * amp + rnorm(1, sd = 0.2)
  }
  fd <- fdata(X, argvals = t)
  list(fd = fd, y = y)
}

# =============================================================================
# fregre.l1
# =============================================================================

test_that("fregre.l1 returns correct structure", {
  d <- .test_robust_data()
  fit <- fregre.l1(d$fd, d$y, ncomp = 3)

  expect_s3_class(fit, "fregre.robust")
  expect_equal(fit$method, "l1")
  expect_equal(fit$ncomp, 3L)
  expect_equal(length(fit$fitted.values), 40)
  expect_equal(length(fit$residuals), 40)
  expect_true(is.numeric(fit$intercept))
  expect_true(is.numeric(fit$beta.t))
  expect_true(is.numeric(fit$coefficients))
  expect_true(is.numeric(fit$r.squared))
  expect_true(fit$r.squared >= 0 && fit$r.squared <= 1)
  expect_true(!is.null(fit$weights))
  expect_true(!is.null(fit$.fpca_mean))
  expect_true(!is.null(fit$.fpca_rotation))
  expect_true(!is.null(fit$.fpca_scores))
})

test_that("fregre.l1 respects ncomp", {
  d <- .test_robust_data()
  fit2 <- fregre.l1(d$fd, d$y, ncomp = 2)
  fit5 <- fregre.l1(d$fd, d$y, ncomp = 5)

  expect_equal(fit2$ncomp, 2L)
  expect_equal(fit5$ncomp, 5L)
  expect_equal(ncol(fit2$.fpca_scores), 2)
  expect_equal(ncol(fit5$.fpca_scores), 5)
})

test_that("fregre.l1 validates input", {
  expect_error(fregre.l1("not_fdata", 1:10), "fdata")
})

# =============================================================================
# fregre.huber
# =============================================================================

test_that("fregre.huber returns correct structure", {
  d <- .test_robust_data()
  fit <- fregre.huber(d$fd, d$y, ncomp = 3)

  expect_s3_class(fit, "fregre.robust")
  expect_equal(fit$method, "huber")
  expect_equal(fit$ncomp, 3L)
  expect_equal(length(fit$fitted.values), 40)
  expect_equal(length(fit$residuals), 40)
  expect_true(is.numeric(fit$r.squared))
  expect_true(!is.null(fit$k))
  expect_equal(fit$k, 1.345)
})

test_that("fregre.huber respects k parameter", {
  d <- .test_robust_data()
  fit_low <- fregre.huber(d$fd, d$y, ncomp = 3, k = 0.5)
  fit_high <- fregre.huber(d$fd, d$y, ncomp = 3, k = 5.0)

  expect_equal(fit_low$k, 0.5)
  expect_equal(fit_high$k, 5.0)
})

test_that("fregre.huber validates input", {
  expect_error(fregre.huber("not_fdata", 1:10), "fdata")
})

# =============================================================================
# predict.fregre.robust
# =============================================================================

test_that("predict.fregre.robust works for l1", {
  d <- .test_robust_data()
  fit <- fregre.l1(d$fd, d$y, ncomp = 3)

  pred <- predict(fit, d$fd[1:5, ])
  expect_equal(length(pred), 5)
  expect_true(is.numeric(pred))
  expect_true(all(is.finite(pred)))
})

test_that("predict.fregre.robust works for huber", {
  d <- .test_robust_data()
  fit <- fregre.huber(d$fd, d$y, ncomp = 3)

  pred <- predict(fit, d$fd[1:5, ])
  expect_equal(length(pred), 5)
  expect_true(is.numeric(pred))
})

test_that("predict.fregre.robust validates input", {
  d <- .test_robust_data()
  fit <- fregre.l1(d$fd, d$y, ncomp = 3)
  expect_error(predict(fit, "not_fdata"), "fdata")
})

# =============================================================================
# print.fregre.robust
# =============================================================================

test_that("print.fregre.robust works for l1", {
  d <- .test_robust_data()
  fit <- fregre.l1(d$fd, d$y, ncomp = 3)
  expect_output(print(fit), "Robust|L1|Huber|fregre", ignore.case = TRUE)
})

test_that("print.fregre.robust works for huber", {
  d <- .test_robust_data()
  fit <- fregre.huber(d$fd, d$y, ncomp = 3)
  expect_output(print(fit), "Robust|L1|Huber|fregre", ignore.case = TRUE)
})

# =============================================================================
# Explainability with robust models
# =============================================================================

test_that("fregre.beta.decomp works with fregre.l1", {
  d <- .test_robust_data()
  fit <- fregre.l1(d$fd, d$y, ncomp = 3)
  result <- fregre.beta.decomp(fit)
  expect_true(!is.null(result))
})

test_that("fregre.pointwise works with fregre.l1", {
  d <- .test_robust_data()
  fit <- fregre.l1(d$fd, d$y, ncomp = 3)
  result <- fregre.pointwise(fit)
  expect_true(!is.null(result))
})

test_that("fregre.shap works with fregre.huber", {
  d <- .test_robust_data()
  fit <- fregre.huber(d$fd, d$y, ncomp = 3)
  result <- fregre.shap(fit, d$fd)
  expect_true(!is.null(result))
  expect_true("values" %in% names(result) || "base_value" %in% names(result))
})

test_that("fregre.influence works with fregre.robust", {
  d <- .test_robust_data()
  fit <- fregre.l1(d$fd, d$y, ncomp = 3)
  result <- fregre.influence(fit, d$fd)
  expect_true(!is.null(result))
  expect_true("leverage" %in% names(result) || "cooks_distance" %in% names(result))
})

test_that("fregre.pdp works with fregre.huber", {
  d <- .test_robust_data()
  fit <- fregre.huber(d$fd, d$y, ncomp = 3)
  result <- fregre.pdp(fit, d$fd, component = 1, n.grid = 10)
  expect_true(!is.null(result))
  expect_true("grid_values" %in% names(result) || "pdp_curve" %in% names(result))
})

test_that("fregre.ale works with fregre.l1", {
  d <- .test_robust_data()
  fit <- fregre.l1(d$fd, d$y, ncomp = 3)
  result <- fregre.ale(fit, d$fd, component = 1, n.bins = 10)
  expect_true(!is.null(result))
  expect_true("ale_values" %in% names(result) || "bin_midpoints" %in% names(result))
})

test_that("fregre.importance works with fregre.robust", {
  d <- .test_robust_data()
  fit <- fregre.huber(d$fd, d$y, ncomp = 3)
  result <- fregre.importance(fit, d$fd, d$y, n.perm = 10, seed = 42)
  expect_true(!is.null(result))
  expect_true("importance" %in% names(result) || "baseline_metric" %in% names(result))
})

test_that("fregre.sobol works with fregre.robust", {
  d <- .test_robust_data()
  fit <- fregre.l1(d$fd, d$y, ncomp = 3)
  result <- fregre.sobol(fit, d$fd, d$y)
  expect_true(!is.null(result))
  expect_true("first_order" %in% names(result) || "total_order" %in% names(result))
})
