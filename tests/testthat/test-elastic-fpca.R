# Tests for elastic FPCA (amplitude, phase, joint)

test_that("vert.fpca returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  X <- matrix(0, 15, 20)
  for (i in 1:15) {
    phase <- runif(1, -0.1, 0.1)
    amp <- rnorm(1, 1, 0.2)
    X[i, ] <- amp * sin(2 * pi * (t + phase))
  }
  fd <- fdata(X, argvals = t)
  km <- karcher.mean(fd, max.iter = 5)

  vfpca <- vert.fpca(km, ncomp = 2)

  expect_s3_class(vfpca, "vert.fpca")
  expect_equal(nrow(vfpca$scores), 15)
  expect_equal(ncol(vfpca$scores), 2)
  expect_equal(vfpca$ncomp, 2)
  expect_equal(length(vfpca$eigenvalues), 2)
  expect_equal(length(vfpca$cumulative.variance), 2)
  expect_true(all(vfpca$eigenvalues >= 0))
  expect_true(all(vfpca$cumulative.variance >= 0))
  expect_true(all(vfpca$cumulative.variance <= 1))
})

test_that("horiz.fpca returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  X <- matrix(0, 15, 20)
  for (i in 1:15) {
    phase <- runif(1, -0.1, 0.1)
    X[i, ] <- sin(2 * pi * (t + phase))
  }
  fd <- fdata(X, argvals = t)
  km <- karcher.mean(fd, max.iter = 5)

  hfpca <- horiz.fpca(km, ncomp = 2)

  expect_s3_class(hfpca, "horiz.fpca")
  expect_equal(nrow(hfpca$scores), 15)
  expect_equal(ncol(hfpca$scores), 2)
  expect_equal(length(hfpca$eigenvalues), 2)
  expect_true(all(hfpca$eigenvalues >= 0))
  expect_true(!is.null(hfpca$eigenfunctions.gam))
})

test_that("joint.fpca returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  X <- matrix(0, 15, 20)
  for (i in 1:15) {
    phase <- runif(1, -0.1, 0.1)
    amp <- rnorm(1, 1, 0.2)
    X[i, ] <- amp * sin(2 * pi * (t + phase))
  }
  fd <- fdata(X, argvals = t)
  km <- karcher.mean(fd, max.iter = 5)

  jfpca <- joint.fpca(km, ncomp = 2)

  expect_s3_class(jfpca, "joint.fpca")
  expect_equal(nrow(jfpca$scores), 15)
  expect_equal(ncol(jfpca$scores), 2)
  expect_true(is.finite(jfpca$balance.c))
  expect_true(jfpca$balance.c > 0)
  expect_true(!is.null(jfpca$vert.component))
  expect_true(!is.null(jfpca$horiz.component))
})

test_that("elastic FPCA cumulative variance is non-decreasing", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  X <- matrix(0, 15, 20)
  for (i in 1:15) {
    phase <- runif(1, -0.1, 0.1)
    amp <- rnorm(1, 1, 0.2)
    X[i, ] <- amp * sin(2 * pi * (t + phase))
  }
  fd <- fdata(X, argvals = t)
  km <- karcher.mean(fd, max.iter = 5)

  vfpca <- vert.fpca(km, ncomp = 3)
  cv <- vfpca$cumulative.variance
  expect_true(all(diff(cv) >= -1e-10))
})

test_that("elastic FPCA input validation works", {
  expect_error(vert.fpca("not_karcher"))
  expect_error(horiz.fpca("not_karcher"))
  expect_error(joint.fpca("not_karcher"))
})

test_that("print methods work for elastic FPCA", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  X <- matrix(0, 15, 20)
  for (i in 1:15) {
    phase <- runif(1, -0.1, 0.1)
    amp <- rnorm(1, 1, 0.2)
    X[i, ] <- amp * sin(2 * pi * (t + phase))
  }
  fd <- fdata(X, argvals = t)
  km <- karcher.mean(fd, max.iter = 5)

  vfpca <- vert.fpca(km, ncomp = 2)
  hfpca <- horiz.fpca(km, ncomp = 2)
  jfpca <- joint.fpca(km, ncomp = 2)

  expect_output(print(vfpca), "Vertical")
  expect_output(print(hfpca), "Horizontal")
  expect_output(print(jfpca), "Joint")
})
