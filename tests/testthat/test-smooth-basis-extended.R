# Extended smooth basis tests: cross-validation, different basis types, prediction

# =============================================================================
# Helper
# =============================================================================

.test_smooth_data <- function(m = 50, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  y_true <- sin(2 * pi * t) + 0.5 * cos(4 * pi * t)
  y_noisy <- y_true + rnorm(m, sd = 0.3)
  fd <- fdata(matrix(y_noisy, 1, m), argvals = t)
  list(fd = fd, y_true = y_true, t = t)
}

# =============================================================================
# smooth.basis.fd with different basis types
# =============================================================================

test_that("smooth.basis.fd with bspline basis", {
  d <- .test_smooth_data()
  result <- smooth.basis.fd(d$fd, type = "bspline", nbasis = 15, lambda = 0.01)

  expect_s3_class(result, "smooth.basis")
  expect_s3_class(result$fitted, "fdata")
  expect_equal(nrow(result$fitted$data), 1)
  expect_true(result$gcv > 0)
})

test_that("smooth.basis.fd with fourier basis", {
  d <- .test_smooth_data()
  result <- smooth.basis.fd(d$fd, type = "fourier", nbasis = 11, lambda = 0.01)

  expect_s3_class(result, "smooth.basis")
  expect_s3_class(result$fitted, "fdata")
})

test_that("smooth.basis.fd with varying lambda", {
  d <- .test_smooth_data()

  fit_low <- smooth.basis.fd(d$fd, nbasis = 15, lambda = 1e-6)
  fit_high <- smooth.basis.fd(d$fd, nbasis = 15, lambda = 100)

  # High lambda should give smoother (lower variance) fit
  range_low <- diff(range(fit_low$fitted$data))
  range_high <- diff(range(fit_high$fitted$data))
  expect_lt(range_high, range_low)
})

test_that("smooth.basis.fd with multiple curves", {
  set.seed(42)
  m <- 30
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, 5, m)
  for (i in 1:5) X[i, ] <- sin(2 * pi * t) + rnorm(m, sd = 0.3)
  fd <- fdata(X, argvals = t)

  result <- smooth.basis.fd(fd, nbasis = 11, lambda = 0.01)
  expect_equal(nrow(result$fitted$data), 5)
  # Residuals can be computed from fitted - original
  resid <- result$fdataobj$data - result$fitted$data
  expect_equal(nrow(resid), 5)
})

test_that("smooth.basis.fd fitted approximates original", {
  d <- .test_smooth_data()
  result <- smooth.basis.fd(d$fd, nbasis = 15, lambda = 0.01)

  # Fitted should be close to original (smoothed version)
  resid <- as.vector(d$fd$data) - as.vector(result$fitted$data)
  expect_true(mean(resid^2) < var(as.vector(d$fd$data)))
})

# =============================================================================
# smooth.basis.gcv
# =============================================================================

test_that("smooth.basis.gcv auto-selects lambda", {
  d <- .test_smooth_data()
  result <- smooth.basis.gcv(d$fd, nbasis = 15)

  expect_s3_class(result, "smooth.basis")
  expect_true(result$lambda > 0)
  expect_true(result$gcv > 0)
})

test_that("smooth.basis.gcv with custom log.lambda.range", {
  d <- .test_smooth_data()
  result <- smooth.basis.gcv(d$fd, nbasis = 15,
                              log.lambda.range = c(-8, 2))

  expect_s3_class(result, "smooth.basis")
  expect_true(result$lambda > 0)
})

test_that("smooth.basis.gcv returns fitted data with correct dimensions", {
  set.seed(42)
  m <- 30
  fd <- fdata(matrix(rnorm(3 * m), 3, m), argvals = seq(0, 1, length.out = m))

  result <- smooth.basis.gcv(fd, nbasis = 11)
  expect_equal(nrow(result$fitted$data), 3)
  expect_equal(ncol(result$fitted$data), m)
})

test_that("smooth.basis.gcv with fourier type", {
  d <- .test_smooth_data()
  result <- smooth.basis.gcv(d$fd, type = "fourier", nbasis = 11)

  expect_s3_class(result, "smooth.basis")
  expect_true(result$lambda > 0)
})

# =============================================================================
# print and plot methods
# =============================================================================

test_that("print.smooth.basis displays info", {
  d <- .test_smooth_data()
  result <- smooth.basis.fd(d$fd, nbasis = 11, lambda = 0.01)

  expect_output(print(result), "Penalized Basis Smoothing")
})

test_that("plot.smooth.basis type='fit' works", {
  d <- .test_smooth_data()
  result <- smooth.basis.fd(d$fd, nbasis = 11, lambda = 0.01)

  pdf(NULL)
  expect_no_error(plot(result, type = "fit"))
  dev.off()
})

test_that("plot.smooth.basis type='residuals' works", {
  d <- .test_smooth_data()
  result <- smooth.basis.fd(d$fd, nbasis = 11, lambda = 0.01)

  pdf(NULL)
  expect_no_error(plot(result, type = "residuals"))
  dev.off()
})
