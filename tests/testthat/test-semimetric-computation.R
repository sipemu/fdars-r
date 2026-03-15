# Tests for semimetric distance computation correctness

# =============================================================================
# Basic distance properties
# =============================================================================

test_that("semimetric.pca gives zero distance for identical curves", {
  set.seed(42)
  X <- matrix(rnorm(50), 5, 10)
  fd <- fdata(X, argvals = seq(0, 1, length.out = 10))
  # Distance of each curve to itself
  D <- semimetric.pca(fd, fd, ncomp = 3)
  expect_true(all(abs(diag(D)) < 1e-10))
})

test_that("semimetric.pca distances are non-negative", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
  D <- semimetric.pca(fd, fd, ncomp = 3)
  expect_true(all(D >= -1e-10))
})

test_that("semimetric.pca self-distance matrix is symmetric", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(50), 5, 10), argvals = seq(0, 1, length.out = 10))
  D <- semimetric.pca(fd, fd, ncomp = 3)
  expect_equal(D, t(D), tolerance = 1e-10)
})

test_that("semimetric.deriv gives zero distance for identical curves", {
  set.seed(42)
  X <- matrix(rnorm(50), 5, 10)
  fd <- fdata(X, argvals = seq(0, 1, length.out = 10))
  D <- semimetric.deriv(fd, fd, nderiv = 1)
  expect_true(all(abs(diag(D)) < 1e-10))
})

test_that("semimetric.deriv distances are non-negative", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
  D <- semimetric.deriv(fd, fd, nderiv = 1)
  expect_true(all(D >= -1e-10))
})

test_that("semimetric.basis gives zero distance for identical curves", {
  set.seed(42)
  X <- matrix(rnorm(50), 5, 10)
  fd <- fdata(X, argvals = seq(0, 1, length.out = 10))
  D <- semimetric.basis(fd, fd, nbasis = 5)
  expect_true(all(abs(diag(D)) < 1e-10))
})

test_that("semimetric.fourier gives zero distance for identical curves", {
  set.seed(42)
  X <- matrix(rnorm(50), 5, 10)
  fd <- fdata(X, argvals = seq(0, 1, length.out = 10))
  D <- semimetric.fourier(fd, fd, nbasis = 5)
  expect_true(all(abs(diag(D)) < 1e-10))
})

test_that("semimetric.hshift gives zero distance for identical curves", {
  set.seed(42)
  X <- matrix(rnorm(50), 5, 10)
  fd <- fdata(X, argvals = seq(0, 1, length.out = 10))
  D <- semimetric.hshift(fd, fd)
  expect_true(all(abs(diag(D)) < 1e-10))
})

# =============================================================================
# Scaled curves should have proportional distances
# =============================================================================

test_that("semimetric.pca scales with amplitude", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  x1 <- sin(2 * pi * t)
  x2 <- sin(2 * pi * t) + 0.5
  x3 <- sin(2 * pi * t) + 2.0

  fd <- fdata(rbind(x1, x2, x3), argvals = t)
  D <- semimetric.pca(fd, fd, ncomp = 2)

  # Distance from x1 to x3 should be larger than x1 to x2
  expect_true(D[1, 3] > D[1, 2])
})

# =============================================================================
# Output dimensions
# =============================================================================

test_that("semimetric.pca returns correct dimensions", {
  set.seed(42)
  fd1 <- fdata(matrix(rnorm(50), 5, 10), argvals = seq(0, 1, length.out = 10))
  fd2 <- fdata(matrix(rnorm(30), 3, 10), argvals = seq(0, 1, length.out = 10))
  D <- semimetric.pca(fd1, fd2, ncomp = 3)
  expect_equal(dim(D), c(5, 3))
})

test_that("semimetric.deriv returns correct dimensions", {
  set.seed(42)
  fd1 <- fdata(matrix(rnorm(80), 8, 10), argvals = seq(0, 1, length.out = 10))
  fd2 <- fdata(matrix(rnorm(40), 4, 10), argvals = seq(0, 1, length.out = 10))
  D <- semimetric.deriv(fd1, fd2, nderiv = 1)
  expect_equal(dim(D), c(8, 4))
})

test_that("semimetric.basis returns correct dimensions", {
  set.seed(42)
  fd1 <- fdata(matrix(rnorm(60), 6, 10), argvals = seq(0, 1, length.out = 10))
  fd2 <- fdata(matrix(rnorm(20), 2, 10), argvals = seq(0, 1, length.out = 10))
  D <- semimetric.basis(fd1, fd2, nbasis = 5)
  expect_equal(dim(D), c(6, 2))
})

test_that("semimetric.fourier returns correct dimensions", {
  set.seed(42)
  fd1 <- fdata(matrix(rnorm(70), 7, 10), argvals = seq(0, 1, length.out = 10))
  fd2 <- fdata(matrix(rnorm(50), 5, 10), argvals = seq(0, 1, length.out = 10))
  D <- semimetric.fourier(fd1, fd2, nbasis = 5)
  expect_equal(dim(D), c(7, 5))
})

test_that("semimetric.hshift returns correct dimensions", {
  set.seed(42)
  fd1 <- fdata(matrix(rnorm(40), 4, 10), argvals = seq(0, 1, length.out = 10))
  fd2 <- fdata(matrix(rnorm(30), 3, 10), argvals = seq(0, 1, length.out = 10))
  D <- semimetric.hshift(fd1, fd2)
  expect_equal(dim(D), c(4, 3))
})
