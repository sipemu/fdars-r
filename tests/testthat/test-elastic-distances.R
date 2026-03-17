# Tests for elastic distance functions, decomposition, constrained alignment, srsf.inverse

# =============================================================================
# Helper
# =============================================================================

.test_elastic_data <- function(n = 5, m = 30, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  for (i in 1:n) {
    phase <- runif(1, -0.1, 0.1)
    X[i, ] <- sin(2 * pi * (t + phase)) + rnorm(m, sd = 0.05)
  }
  fdata(X, argvals = t)
}

# =============================================================================
# amplitude.distance
# =============================================================================

test_that("amplitude.distance returns symmetric distance matrix", {
  fd <- .test_elastic_data(n = 4, m = 20)
  D <- amplitude.distance(fd)

  expect_true(is.matrix(D))
  expect_equal(nrow(D), 4)
  expect_equal(ncol(D), 4)
  expect_equal(D, t(D), tolerance = 1e-6)
  expect_true(all(diag(D) < 1e-6))
  expect_true(all(D >= -1e-10))
})

test_that("amplitude.distance identical curves have zero distance", {
  t <- seq(0, 1, length.out = 20)
  X <- matrix(sin(2 * pi * t), 3, 20, byrow = TRUE)
  fd <- fdata(X, argvals = t)

  D <- amplitude.distance(fd)
  expect_true(all(D < 1e-6))
})

test_that("amplitude.distance with lambda parameter", {
  fd <- .test_elastic_data(n = 4, m = 20)
  D0 <- amplitude.distance(fd, lambda = 0)
  D1 <- amplitude.distance(fd, lambda = 1)

  expect_true(is.matrix(D0))
  expect_true(is.matrix(D1))
})

test_that("amplitude.distance with reference data", {
  fd <- .test_elastic_data(n = 4, m = 20)
  fd_ref <- .test_elastic_data(n = 3, m = 20, seed = 99)

  D <- amplitude.distance(fd, fdataref = fd_ref)
  expect_equal(nrow(D), 4)
  expect_equal(ncol(D), 3)
})

test_that("amplitude.distance validates input", {
  expect_error(amplitude.distance("not_fdata"), "fdata")
})

# =============================================================================
# phase.distance
# =============================================================================

test_that("phase.distance returns symmetric distance matrix", {
  fd <- .test_elastic_data(n = 4, m = 20)
  D <- phase.distance(fd)

  expect_true(is.matrix(D))
  expect_equal(nrow(D), 4)
  expect_equal(ncol(D), 4)
  expect_equal(D, t(D), tolerance = 1e-6)
  expect_true(all(diag(D) < 1e-6))
  expect_true(all(D >= -1e-10))
})

test_that("phase.distance identical curves have zero distance", {
  t <- seq(0, 1, length.out = 20)
  X <- matrix(sin(2 * pi * t), 3, 20, byrow = TRUE)
  fd <- fdata(X, argvals = t)

  D <- phase.distance(fd)
  expect_true(all(D < 1e-6))
})

test_that("phase.distance validates input", {
  expect_error(phase.distance("not_fdata"), "fdata")
})

# =============================================================================
# elastic.decomposition
# =============================================================================

test_that("elastic.decomposition returns amplitude and phase components", {
  t <- seq(0, 1, length.out = 30)
  f1 <- sin(2 * pi * t)
  f2 <- sin(2 * pi * (t + 0.05))

  result <- elastic.decomposition(f1, f2, argvals = t)

  expect_true(is.list(result))
  expect_true("d_amplitude" %in% names(result))
  expect_true("d_phase" %in% names(result))
  expect_true(is.numeric(result$d_amplitude))
  expect_true(is.numeric(result$d_phase))
  expect_gte(result$d_amplitude, 0)
  expect_gte(result$d_phase, 0)
})

test_that("elastic.decomposition with fdata input", {
  t <- seq(0, 1, length.out = 30)
  f1 <- fdata(matrix(sin(2 * pi * t), 1, 30), argvals = t)
  f2 <- fdata(matrix(cos(2 * pi * t), 1, 30), argvals = t)

  result <- elastic.decomposition(f1, f2)
  expect_true(is.list(result))
  expect_gte(result$d_amplitude, 0)
  expect_gte(result$d_phase, 0)
})

test_that("elastic.decomposition identical curves have zero components", {
  t <- seq(0, 1, length.out = 30)
  f <- sin(2 * pi * t)

  result <- elastic.decomposition(f, f, argvals = t)
  expect_lt(result$d_amplitude, 1e-4)
  expect_lt(result$d_phase, 1e-4)
})

test_that("elastic.decomposition validates argvals requirement", {
  expect_error(elastic.decomposition(1:10, 1:10), "argvals")
})

# =============================================================================
# elastic.align.constrained
# =============================================================================

test_that("elastic.align.constrained returns aligned data and warps", {
  fd <- .test_elastic_data(n = 5, m = 30)
  result <- elastic.align.constrained(fd)

  expect_s3_class(result, "elastic.align")
  expect_s3_class(result$aligned, "fdata")
  expect_s3_class(result$gammas, "fdata")
  expect_equal(nrow(result$aligned$data), 5)
  expect_equal(length(result$distances), 5)
  expect_true(all(result$distances >= 0))
})

test_that("elastic.align.constrained with explicit target", {
  t <- seq(0, 1, length.out = 30)
  fd <- .test_elastic_data(n = 5, m = 30)
  target <- fdata(matrix(sin(2 * pi * t), 1, 30), argvals = t)

  result <- elastic.align.constrained(fd, target = target)
  expect_true(is.list(result))
})

test_that("elastic.align.constrained validates input", {
  expect_error(elastic.align.constrained("not_fdata"), "fdata")
})

# =============================================================================
# srsf.inverse
# =============================================================================

test_that("srsf.inverse returns numeric vector", {
  t <- seq(0, 1, length.out = 30)
  q <- cos(2 * pi * t)  # an SRSF

  f <- srsf.inverse(q, t)
  expect_true(is.numeric(f))
  expect_equal(length(f), 30)
  expect_true(all(is.finite(f)))
})

test_that("srsf.inverse with custom f0", {
  t <- seq(0, 1, length.out = 30)
  q <- rep(1, 30)  # constant SRSF

  f1 <- srsf.inverse(q, t, f0 = 0)
  f2 <- srsf.inverse(q, t, f0 = 5)

  # Different starting values should give different results
  expect_false(isTRUE(all.equal(f1, f2)))
  # But shifted by f0 difference
  expect_equal(f2 - f1, rep(5, 30), tolerance = 1e-6)
})

test_that("srsf roundtrip: srsf.transform -> inverse -> srsf.transform", {
  t <- seq(0, 1, length.out = 50)
  f_orig <- sin(2 * pi * t) + 1

  # Forward: f -> q
  fd <- fdata(matrix(f_orig, 1, 50), argvals = t)
  q_fd <- srsf.transform(fd)
  q <- as.numeric(q_fd$data)

  # Inverse: q -> f
  f_recovered <- srsf.inverse(q, t, f0 = f_orig[1])

  # Should approximately recover original
  expect_equal(f_recovered, f_orig, tolerance = 0.1)
})

# =============================================================================
# elastic.attribution
# =============================================================================

test_that("elastic.attribution returns attribution results", {
  set.seed(42)
  m <- 20
  t <- seq(0, 1, length.out = m)
  n <- 30
  X <- matrix(0, n, m)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(m, sd = 0.1)
    y[i] <- amp + rnorm(1, sd = 0.1)
  }
  fd <- fdata(X, argvals = t)

  result <- elastic.attribution(fd, y, ncomp = 2, n.perm = 10, seed = 42)

  expect_true(is.list(result))
})
