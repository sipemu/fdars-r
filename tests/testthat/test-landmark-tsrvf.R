# Tests for landmark registration and TSRVF functions

# =============================================================================
# landmark.register
# =============================================================================

test_that("landmark.register returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 50)
  n <- 10
  X <- matrix(0, n, 50)
  for (i in 1:n) {
    shift <- rnorm(1, 0, 0.05)
    X[i, ] <- sin(2 * pi * (t - shift))
  }
  fd <- fdata(X, argvals = t)

  result <- landmark.register(fd, kind = "peak", min.prominence = 0.3)

  expect_s3_class(result, "landmark.register")
  expect_s3_class(result$registered, "fdata")
  expect_s3_class(result$gammas, "fdata")
  expect_equal(nrow(result$registered$data), n)
  expect_equal(nrow(result$gammas$data), n)
  expect_true(!is.null(result$target_landmarks))
})

test_that("landmark.register works with valley kind", {
  set.seed(42)
  t <- seq(0, 1, length.out = 50)
  n <- 8
  X <- matrix(0, n, 50)
  for (i in 1:n) {
    shift <- rnorm(1, 0, 0.05)
    X[i, ] <- sin(2 * pi * (t - shift))
  }
  fd <- fdata(X, argvals = t)

  result <- landmark.register(fd, kind = "valley", min.prominence = 0.3)

  expect_s3_class(result, "landmark.register")
})

test_that("landmark.register validates input", {
  expect_error(landmark.register("not_fdata"), "fdata")
})

test_that("print.landmark.register works", {
  set.seed(42)
  t <- seq(0, 1, length.out = 50)
  X <- matrix(0, 8, 50)
  for (i in 1:8) X[i, ] <- sin(2 * pi * (t - rnorm(1, 0, 0.05)))
  fd <- fdata(X, argvals = t)

  result <- landmark.register(fd, kind = "peak", min.prominence = 0.3)
  expect_output(print(result), "Landmark Registration")
})

test_that("plot.landmark.register runs without error", {
  set.seed(42)
  t <- seq(0, 1, length.out = 50)
  X <- matrix(0, 8, 50)
  for (i in 1:8) X[i, ] <- sin(2 * pi * (t - rnorm(1, 0, 0.05)))
  fd <- fdata(X, argvals = t)

  result <- landmark.register(fd, kind = "peak", min.prominence = 0.3)

  pdf(NULL)
  expect_no_error(plot(result, type = "registered"))
  expect_no_error(plot(result, type = "original"))
  expect_no_error(plot(result, type = "warps"))
  dev.off()
})

# =============================================================================
# tsrvf.from.alignment / tsrvf.inverse
# =============================================================================

test_that("tsrvf.from.alignment returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 30)
  n <- 15
  X <- matrix(0, n, 30)
  for (i in 1:n) {
    phase <- rnorm(1, 0, 0.05)
    X[i, ] <- sin(2 * pi * (t + phase)) + rnorm(30, sd = 0.05)
  }
  fd <- fdata(X, argvals = t)
  km <- karcher.mean(fd, max.iter = 3)

  result <- tsrvf.from.alignment(km)

  expect_s3_class(result, "tsrvf")
  expect_s3_class(result$tangent_vectors, "fdata")
  expect_equal(nrow(result$tangent_vectors$data), n)
  expect_true(is.numeric(result$mean_srsf))
  expect_true(is.numeric(result$mean_srsf_norm))
  expect_true(is.logical(result$converged) || is.numeric(result$converged))
})

test_that("tsrvf.inverse reconstructs curves", {
  set.seed(42)
  t <- seq(0, 1, length.out = 30)
  n <- 10
  X <- matrix(0, n, 30)
  for (i in 1:n) {
    phase <- rnorm(1, 0, 0.03)
    X[i, ] <- sin(2 * pi * (t + phase)) + rnorm(30, sd = 0.03)
  }
  fd <- fdata(X, argvals = t)
  km <- karcher.mean(fd, max.iter = 3)
  tv <- tsrvf.from.alignment(km)

  recon <- tsrvf.inverse(tv)

  expect_s3_class(recon, "fdata")
  expect_equal(nrow(recon$data), n)
  expect_equal(ncol(recon$data), 30)
})

test_that("tsrvf.inverse validates input", {
  expect_error(tsrvf.inverse("not_tsrvf"), "tsrvf")
})

test_that("print.tsrvf works", {
  set.seed(42)
  t <- seq(0, 1, length.out = 30)
  X <- matrix(0, 10, 30)
  for (i in 1:10) X[i, ] <- sin(2 * pi * t) + rnorm(30, sd = 0.05)
  fd <- fdata(X, argvals = t)
  km <- karcher.mean(fd, max.iter = 3)
  tv <- tsrvf.from.alignment(km)

  expect_output(print(tv), "TSRVF")
})

test_that("plot.tsrvf runs for all types", {
  set.seed(42)
  t <- seq(0, 1, length.out = 30)
  X <- matrix(0, 10, 30)
  for (i in 1:10) X[i, ] <- sin(2 * pi * t) + rnorm(30, sd = 0.05)
  fd <- fdata(X, argvals = t)
  km <- karcher.mean(fd, max.iter = 3)
  tv <- tsrvf.from.alignment(km)

  pdf(NULL)
  expect_no_error(plot(tv, type = "tangent"))
  expect_no_error(plot(tv, type = "mean"))
  expect_no_error(plot(tv, type = "gammas"))
  dev.off()
})
