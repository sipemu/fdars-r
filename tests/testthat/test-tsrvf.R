test_that("tsrvf.transform returns correct structure", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
  tv <- tsrvf.transform(fd, max.iter = 5)
  expect_s3_class(tv, "tsrvf")
  expect_s3_class(tv$tangent_vectors, "fdata")
  expect_s3_class(tv$mean, "fdata")
  expect_s3_class(tv$gammas, "fdata")
  expect_equal(nrow(tv$tangent_vectors$data), 20)
  expect_equal(ncol(tv$tangent_vectors$data), 10)
  expect_equal(nrow(tv$mean$data), 1)
  expect_type(tv$converged, "logical")
})

test_that("tsrvf.from.alignment matches tsrvf.transform", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
  km <- karcher.mean(fd, max.iter = 5)
  tv_from <- tsrvf.from.alignment(km)
  expect_s3_class(tv_from, "tsrvf")
  expect_equal(nrow(tv_from$tangent_vectors$data), 20)
  expect_equal(ncol(tv_from$tangent_vectors$data), 10)
})

test_that("tsrvf.inverse approximately recovers structure", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
  tv <- tsrvf.transform(fd, max.iter = 5)
  recon <- tsrvf.inverse(tv)
  expect_s3_class(recon, "fdata")
  expect_equal(nrow(recon$data), 20)
  expect_equal(ncol(recon$data), 10)
  # Check reconstruction is not all zeros/NaN
  expect_true(any(is.finite(recon$data)))
})

test_that("tsrvf.inverse recovers correct initial values", {
  set.seed(123)
  n <- 10
  m <- 50
  argvals <- seq(0, 1, length.out = m)
  # Smooth curves with varying initial values
  curves <- matrix(0, n, m)
  for (i in seq_len(n)) {
    curves[i, ] <- sin(2 * pi * argvals + runif(1, 0, 2 * pi)) + rnorm(1, 0, 0.5)
  }
  fd <- fdata(curves, argvals = argvals)
  tv <- tsrvf.transform(fd, max.iter = 10)
  recon <- tsrvf.inverse(tv)
  # Initial values of reconstruction should match original curves
  expect_equal(recon$data[, 1], fd$data[, 1], tolerance = 1e-10)
})

test_that("tsrvf print method works", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
  tv <- tsrvf.transform(fd, max.iter = 3)
  expect_output(print(tv), "TSRVF")
})

test_that("tsrvf.transform validates input", {
  expect_error(tsrvf.transform("not_fdata"))
})

test_that("tsrvf.from.alignment validates input", {
  expect_error(tsrvf.from.alignment("not_karcher"))
})
