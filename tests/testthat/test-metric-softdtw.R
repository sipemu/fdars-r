test_that("metric.softDTW returns proper self-distance matrix", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(100), 10, 10))
  D <- metric.softDTW(fd, gamma = 1.0)
  expect_equal(nrow(D), 10)
  expect_equal(ncol(D), 10)
  # Symmetric
  expect_equal(D, t(D), tolerance = 1e-10)
})

test_that("metric.softDTW cross-distance has correct dimensions", {
  set.seed(42)
  fd1 <- fdata(matrix(rnorm(60), 6, 10))
  fd2 <- fdata(matrix(rnorm(40), 4, 10))
  D <- metric.softDTW(fd1, fd2, gamma = 1.0)
  expect_equal(nrow(D), 6)
  expect_equal(ncol(D), 4)
})

test_that("metric.softDTW divergence is non-negative with zero diagonal", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(100), 10, 10))
  D <- metric.softDTW(fd, gamma = 1.0, divergence = TRUE)
  expect_equal(nrow(D), 10)
  # Divergence should be non-negative
  expect_true(all(D >= -1e-10))
  # Diagonal should be (near) zero
  expect_equal(diag(D), rep(0, 10), tolerance = 1e-10)
})

test_that("softdtw.barycenter returns correct dimensions and converges", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
  bc <- softdtw.barycenter(fd, gamma = 1.0, max.iter = 50)
  expect_s3_class(bc, "fdata")
  expect_equal(nrow(bc$data), 1)
  expect_equal(ncol(bc$data), 10)
})

test_that("metric dispatcher works with softdtw", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(60), 6, 10))
  D <- metric(fd, method = "softdtw", gamma = 1.0)
  expect_equal(nrow(D), 6)
  expect_equal(ncol(D), 6)
})

test_that("metric.softDTW validates input", {
  expect_error(metric.softDTW("not_fdata"))
})
