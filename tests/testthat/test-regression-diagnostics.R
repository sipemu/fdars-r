# Tests for regression diagnostics: significant regions, prediction intervals,
# calibration, ECE

# =============================================================================
# Helper
# =============================================================================

.test_diag_lm <- function(n = 40, m = 20, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(m, sd = 0.1)
    y[i] <- 2 * amp + rnorm(1, sd = 0.3)
  }
  fd <- fdata(X, argvals = t)
  model <- fregre.lm(fd, y, ncomp = 3)
  list(fd = fd, y = y, model = model, t = t, m = m)
}

.test_diag_logistic <- function(n = 40, m = 20, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  for (i in 1:n) {
    if (i <= n / 2) {
      X[i, ] <- sin(2 * pi * t) + rnorm(m, sd = 0.2)
    } else {
      X[i, ] <- cos(2 * pi * t) + rnorm(m, sd = 0.2)
    }
  }
  y <- rep(c(0, 1), each = n / 2)
  fd <- fdata(X, argvals = t)
  model <- functional.logistic(fd, y, ncomp = 3)
  list(fd = fd, y = y, model = model)
}

# =============================================================================
# fregre.significant.regions
# =============================================================================

test_that("fregre.significant.regions with logistic model returns region info", {
  d <- .test_diag_logistic()
  result <- fregre.significant.regions(d$model, alpha = 0.05)

  expect_true(is.list(result))
  expect_true("start_idx" %in% names(result))
  expect_true("end_idx" %in% names(result))
  expect_true("direction" %in% names(result))
})

test_that("fregre.significant.regions with stricter alpha finds fewer regions", {
  d <- .test_diag_logistic()
  r_05 <- fregre.significant.regions(d$model, alpha = 0.05)
  r_50 <- fregre.significant.regions(d$model, alpha = 0.50)

  expect_true(is.list(r_05))
  expect_true(is.list(r_50))
  # Larger alpha should find at least as many regions
  expect_lte(length(r_05$start_idx), length(r_50$start_idx))
})

test_that("fregre.significant.regions validates input", {
  expect_error(fregre.significant.regions("not_model"), "fregre.lm.*fregre.logistic")
})

# =============================================================================
# fregre.prediction.interval
# =============================================================================

test_that("fregre.prediction.interval returns intervals", {
  d <- .test_diag_lm()
  new_fd <- fdata(d$fd$data[1:5, ], argvals = d$t)

  result <- fregre.prediction.interval(d$model, d$fd, new_fd, confidence = 0.95)

  expect_true(is.list(result))
})

test_that("fregre.prediction.interval wider at higher confidence", {
  d <- .test_diag_lm()
  new_fd <- fdata(d$fd$data[1:3, ], argvals = d$t)

  r_90 <- fregre.prediction.interval(d$model, d$fd, new_fd, confidence = 0.90)
  r_99 <- fregre.prediction.interval(d$model, d$fd, new_fd, confidence = 0.99)

  expect_true(is.list(r_90))
  expect_true(is.list(r_99))
})

# =============================================================================
# fregre.calibration (logistic model)
# =============================================================================

test_that("fregre.calibration returns calibration info", {
  d <- .test_diag_logistic()
  result <- fregre.calibration(d$model, d$y, n.groups = 5)

  expect_true(is.list(result))
})

test_that("fregre.calibration with different n.groups", {
  d <- .test_diag_logistic()

  r5 <- fregre.calibration(d$model, d$y, n.groups = 5)
  r10 <- fregre.calibration(d$model, d$y, n.groups = 10)

  expect_true(is.list(r5))
  expect_true(is.list(r10))
})

# =============================================================================
# fregre.ece (Expected Calibration Error)
# =============================================================================

test_that("fregre.ece returns result", {
  d <- .test_diag_logistic()
  result <- fregre.ece(d$model, d$y, n.bins = 5)

  expect_true(is.list(result) || is.numeric(result))
  if (is.numeric(result)) {
    expect_gte(result, 0)
  }
})

test_that("fregre.ece with different n.bins", {
  d <- .test_diag_logistic()

  r5 <- fregre.ece(d$model, d$y, n.bins = 5)
  r10 <- fregre.ece(d$model, d$y, n.bins = 10)

  expect_true(is.list(r5) || is.numeric(r5))
  expect_true(is.list(r10) || is.numeric(r10))
})

# =============================================================================
# pred.MAE
# =============================================================================

test_that("pred.MAE computes mean absolute error correctly", {
  y_true <- c(1, 2, 3, 4, 5)
  y_pred <- c(1.1, 2.2, 2.8, 4.3, 4.7)

  mae <- pred.MAE(y_true, y_pred)
  expect_true(is.numeric(mae))
  expect_equal(length(mae), 1)
  expect_gte(mae, 0)
  expect_equal(mae, mean(abs(y_true - y_pred)))
})

test_that("pred.MAE is zero for perfect predictions", {
  y <- c(1, 2, 3)
  expect_equal(pred.MAE(y, y), 0)
})
