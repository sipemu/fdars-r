test_that("srsf.transform returns correct dimensions", {
  fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
  q <- srsf.transform(fd)
  expect_s3_class(q, "fdata")
  expect_equal(nrow(q$data), 20)
  expect_equal(ncol(q$data), 10)
})

test_that("srsf round-trip approximately recovers original", {
  set.seed(42)
  argvals <- seq(0, 1, length.out = 50)
  f <- sin(2 * pi * argvals)
  fd <- fdata(matrix(f, nrow = 1), argvals = argvals)
  q <- srsf.transform(fd)
  f_hat <- srsf.inverse(q$data[1, ], argvals, f[1])
  # Approximate recovery (integration introduces small error)
  expect_true(cor(f, f_hat) > 0.95)
})

test_that("elastic.align returns correct structure", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(100), 10, 10), argvals = seq(0, 1, length.out = 10))
  res <- elastic.align(fd)
  expect_s3_class(res, "elastic.align")
  expect_s3_class(res$aligned, "fdata")
  expect_s3_class(res$gammas, "fdata")
  expect_equal(nrow(res$aligned$data), 10)
  expect_equal(length(res$distances), 10)
})

test_that("elastic.distance returns symmetric matrix with zero diagonal", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(60), 6, 10), argvals = seq(0, 1, length.out = 10))
  D <- elastic.distance(fd)
  expect_equal(nrow(D), 6)
  expect_equal(ncol(D), 6)
  # Symmetric
  expect_equal(D, t(D), tolerance = 1e-10)
  # Zero diagonal
  expect_equal(diag(D), rep(0, 6), tolerance = 1e-10)
})

test_that("elastic.distance cross-distance has correct dimensions", {
  set.seed(42)
  fd1 <- fdata(matrix(rnorm(60), 6, 10), argvals = seq(0, 1, length.out = 10))
  fd2 <- fdata(matrix(rnorm(40), 4, 10), argvals = seq(0, 1, length.out = 10))
  D <- elastic.distance(fd1, fd2)
  expect_equal(nrow(D), 6)
  expect_equal(ncol(D), 4)
})

test_that("karcher.mean returns correct structure", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(100), 10, 10), argvals = seq(0, 1, length.out = 10))
  km <- karcher.mean(fd, max.iter = 5)
  expect_s3_class(km, "karcher.mean")
  expect_s3_class(km$mean, "fdata")
  expect_equal(nrow(km$mean$data), 1)
  expect_equal(ncol(km$mean$data), 10)
  expect_true(km$n.iter >= 1)
  expect_type(km$converged, "logical")
})

test_that("metric dispatcher works with elastic", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(60), 6, 10), argvals = seq(0, 1, length.out = 10))
  D <- metric(fd, method = "elastic")
  expect_equal(nrow(D), 6)
  expect_equal(ncol(D), 6)
})

test_that("alignment input validation works", {
  expect_error(srsf.transform("not_fdata"))
  expect_error(elastic.align("not_fdata"))
  expect_error(elastic.distance("not_fdata"))
  expect_error(karcher.mean("not_fdata"))
})

test_that("print methods work without error", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(100), 10, 10), argvals = seq(0, 1, length.out = 10))
  res <- elastic.align(fd)
  expect_output(print(res), "Elastic Alignment")

  km <- karcher.mean(fd, max.iter = 3)
  expect_output(print(km), "Karcher Mean")
})
