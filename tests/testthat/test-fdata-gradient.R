# Tests for fdata.gradient and related fdata operations

# =============================================================================
# fdata.gradient
# =============================================================================

test_that("fdata.gradient returns fdata object", {
  t <- seq(0, 1, length.out = 50)
  X <- matrix(sin(2 * pi * t), nrow = 1)
  fd <- fdata(X, argvals = t)

  gd <- fdata.gradient(fd)

  expect_s3_class(gd, "fdata")
  expect_equal(nrow(gd$data), 1)
  expect_equal(ncol(gd$data), 50)
  expect_equal(as.numeric(gd$argvals), t)
})

test_that("fdata.gradient approximates analytical derivative", {
  t <- seq(0, 1, length.out = 100)
  X <- matrix(sin(2 * pi * t), nrow = 1)
  fd <- fdata(X, argvals = t)

  gd <- fdata.gradient(fd)
  expected <- 2 * pi * cos(2 * pi * t)

  # Numerical gradient should be close to analytical (except at boundaries)
  mid <- 10:90
  correlation <- cor(as.numeric(gd$data[1, mid]), expected[mid])
  expect_true(correlation > 0.99)
})

test_that("fdata.gradient handles multiple curves", {
  set.seed(42)
  t <- seq(0, 1, length.out = 30)
  X <- matrix(rnorm(5 * 30), 5, 30)
  fd <- fdata(X, argvals = t)

  gd <- fdata.gradient(fd)

  expect_equal(nrow(gd$data), 5)
  expect_equal(ncol(gd$data), 30)
})

test_that("fdata.gradient preserves names", {
  t <- seq(0, 1, length.out = 20)
  fd <- fdata(matrix(rnorm(20), 1, 20), argvals = t,
              names = list(main = "test", xlab = "time", ylab = "value"))

  gd <- fdata.gradient(fd)

  expect_equal(gd$names$xlab, "time")
  expect_true(grepl("D\\[", gd$names$ylab))
})

test_that("fdata.gradient of constant is near zero", {
  t <- seq(0, 1, length.out = 30)
  X <- matrix(5, nrow = 1, ncol = 30)
  fd <- fdata(X, argvals = t)

  gd <- fdata.gradient(fd)

  expect_true(max(abs(gd$data)) < 1e-10)
})

test_that("fdata.gradient validates input", {
  expect_error(fdata.gradient("not_fdata"), "fdata")
})

# =============================================================================
# model.selection.ncomp
# =============================================================================

test_that("model.selection.ncomp returns optimal components", {
  set.seed(42)
  t <- seq(0, 1, length.out = 15)
  n <- 40
  X <- matrix(0, n, 15)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1, 1, 0.3)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(15, sd = 0.1)
    y[i] <- 2 * amp + rnorm(1, sd = 0.3)
  }
  fd <- fdata(X, argvals = t)

  result <- model.selection.ncomp(fd, y, max.ncomp = 5, criterion = "bic")

  expect_true(!is.null(result$best.ncomp))
  expect_true(result$best.ncomp >= 1 && result$best.ncomp <= 5)
  expect_true(!is.null(result$bic))
  expect_true(!is.null(result$aic))
  expect_true(!is.null(result$gcv))
})

test_that("model.selection.ncomp works with all criteria", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(400), 40, 10), argvals = seq(0, 1, length.out = 10))
  y <- rnorm(40)

  r_aic <- model.selection.ncomp(fd, y, max.ncomp = 5, criterion = "aic")
  r_bic <- model.selection.ncomp(fd, y, max.ncomp = 5, criterion = "bic")
  r_gcv <- model.selection.ncomp(fd, y, max.ncomp = 5, criterion = "gcv")

  expect_true(r_aic$best.ncomp >= 1)
  expect_true(r_bic$best.ncomp >= 1)
  expect_true(r_gcv$best.ncomp >= 1)
})

test_that("model.selection.ncomp validates input", {
  expect_error(model.selection.ncomp("not_fdata", 1:10), "fdata")
})

# =============================================================================
# basis2fdata_2d
# =============================================================================

test_that("fdata2basis_2d and basis2fdata_2d round-trip", {
  set.seed(42)
  m_s <- 10; m_t <- 8
  s <- seq(0, 1, length.out = m_s)
  t <- seq(0, 1, length.out = m_t)
  surface <- outer(sin(2 * pi * s), cos(2 * pi * t))
  fd2d <- fdata(array(surface, dim = c(1, m_s, m_t)))

  coefs <- fdata2basis_2d(fd2d, nbasis.s = 5, nbasis.t = 5, type = "bspline")
  expect_true(is.matrix(coefs))
  expect_equal(nrow(coefs), 1)

  argvals <- list(s, t)
  reconstructed <- basis2fdata_2d(coefs, argvals, nbasis.s = 5, nbasis.t = 5,
                                   type = "bspline")

  expect_true(inherits(reconstructed, "fdata"))
  expect_equal(nrow(reconstructed$data), 1)
})
