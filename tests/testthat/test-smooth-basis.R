# Tests for penalized basis smoothing

test_that("smooth.basis.fd returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 50)
  X <- matrix(0, 10, 50)
  for (i in 1:10) X[i, ] <- sin(2 * pi * t) + rnorm(50, sd = 0.2)
  fd <- fdata(X, argvals = t)

  result <- smooth.basis.fd(fd, type = "bspline", nbasis = 15, lambda = 0.01)

  expect_s3_class(result, "smooth.basis")
  expect_s3_class(result$fitted, "fdata")
  expect_equal(nrow(result$fitted$data), 10)
  expect_equal(ncol(result$fitted$data), 50)
  expect_equal(result$type, "bspline")
  expect_equal(result$nbasis, 15)
  expect_true(is.finite(result$gcv))
  expect_true(is.finite(result$edf))
  expect_true(is.finite(result$aic))
  expect_true(is.finite(result$bic))
  expect_true(!is.null(result$coefficients))
})

test_that("smooth.basis.fd works with fourier basis", {
  set.seed(42)
  t <- seq(0, 1, length.out = 50)
  X <- matrix(sin(2 * pi * t) + rnorm(50, sd = 0.1), nrow = 1)
  fd <- fdata(X, argvals = t)

  result <- smooth.basis.fd(fd, type = "fourier", nbasis = 11, lambda = 0.001)

  expect_s3_class(result, "smooth.basis")
  expect_equal(result$type, "fourier")
  expect_true(is.finite(result$gcv))
})

test_that("smooth.basis.fd smoothing reduces noise", {
  set.seed(42)
  t <- seq(0, 1, length.out = 100)
  signal <- sin(2 * pi * t)
  noisy <- signal + rnorm(100, sd = 0.3)
  fd <- fdata(matrix(noisy, nrow = 1), argvals = t)

  result <- smooth.basis.fd(fd, type = "fourier", nbasis = 11, lambda = 0.001)

  # Smoothed should be closer to true signal than noisy data
  resid_noisy <- sum((noisy - signal)^2)
  resid_smooth <- sum((result$fitted$data[1, ] - signal)^2)
  expect_lt(resid_smooth, resid_noisy)
})

test_that("smooth.basis.gcv selects lambda automatically", {
  set.seed(42)
  t <- seq(0, 1, length.out = 50)
  X <- matrix(0, 5, 50)
  for (i in 1:5) X[i, ] <- sin(4 * pi * t) + rnorm(50, sd = 0.2)
  fd <- fdata(X, argvals = t)

  result <- smooth.basis.gcv(fd, type = "bspline", nbasis = 15)

  expect_s3_class(result, "smooth.basis")
  expect_s3_class(result$fitted, "fdata")
  expect_true(is.finite(result$gcv))
  expect_true(is.finite(result$edf))
})

test_that("smooth.basis.fd input validation works", {
  expect_error(smooth.basis.fd("not_fdata"))
})

test_that("print.smooth.basis works without error", {
  set.seed(42)
  t <- seq(0, 1, length.out = 50)
  fd <- fdata(matrix(sin(2 * pi * t) + rnorm(50, sd = 0.1), nrow = 1), argvals = t)
  result <- smooth.basis.fd(fd, type = "bspline", nbasis = 10, lambda = 0.01)

  expect_output(print(result), "Penalized Basis Smoothing")
  expect_output(print(result), "bspline")
})

test_that("plot.smooth.basis works without error", {
  set.seed(42)
  t <- seq(0, 1, length.out = 50)
  X <- matrix(0, 5, 50)
  for (i in 1:5) X[i, ] <- sin(2 * pi * t) + rnorm(50, sd = 0.1)
  fd <- fdata(X, argvals = t)
  result <- smooth.basis.fd(fd, type = "bspline", nbasis = 10, lambda = 0.01)

  expect_no_error(plot(result, type = "fit"))
  expect_no_error(plot(result, type = "residuals"))
})
