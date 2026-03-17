# Tests for Phase 8: Shape analysis functions

# =============================================================================
# shape.representative
# =============================================================================

test_that("shape.representative returns correct structure", {
  t <- seq(0, 1, length.out = 30)
  f <- sin(2 * pi * t)

  result <- shape.representative(f, argvals = t)

  expect_true(is.list(result))
  expect_true(is.numeric(result$representative))
  expect_equal(length(result$representative), 30)
  expect_true(is.numeric(result$representative.srsf))
  expect_equal(length(result$representative.srsf), 30)
  expect_true(is.numeric(result$gamma))
  expect_equal(length(result$gamma), 30)
  expect_true(is.numeric(result$translation))
  expect_true(is.numeric(result$scale))
})

test_that("shape.representative reparameterization quotient", {
  t <- seq(0, 1, length.out = 30)
  f <- sin(2 * pi * t)

  result <- shape.representative(f, argvals = t,
                                   quotient = "reparameterization")

  expect_true(is.numeric(result$representative))
  # For reparameterization, translation should be 0 and scale should be 1
  expect_equal(result$translation, 0, tolerance = 1e-6)
  expect_equal(result$scale, 1, tolerance = 1e-6)
})

test_that("shape.representative translation quotient", {
  t <- seq(0, 1, length.out = 30)
  f <- sin(2 * pi * t) + 5

  result <- shape.representative(f, argvals = t, quotient = "translation")

  expect_true(is.numeric(result$representative))
  expect_true(is.numeric(result$translation))
})

test_that("shape.representative scale quotient", {
  t <- seq(0, 1, length.out = 30)
  f <- 3 * sin(2 * pi * t)

  result <- shape.representative(f, argvals = t, quotient = "scale")

  expect_true(is.numeric(result$representative))
  expect_true(is.numeric(result$scale))
})

test_that("shape.representative gamma is valid warp", {
  t <- seq(0, 1, length.out = 30)
  f <- sin(2 * pi * t)

  result <- shape.representative(f, argvals = t)

  # Boundary conditions
  expect_equal(result$gamma[1], 0, tolerance = 1e-4)
  expect_equal(result$gamma[30], 1, tolerance = 1e-4)
  # Monotonicity
  diffs <- diff(result$gamma)
  expect_true(all(diffs >= -1e-10))
})

test_that("shape.representative default argvals", {
  f <- sin(2 * pi * seq(0, 1, length.out = 20))

  result <- shape.representative(f)

  expect_true(is.list(result))
  expect_equal(length(result$representative), 20)
})

# =============================================================================
# shape.distance
# =============================================================================

test_that("shape.distance between identical curves is near zero", {
  t <- seq(0, 1, length.out = 30)
  f <- sin(2 * pi * t)

  result <- shape.distance(f, f, argvals = t)

  expect_true(is.list(result))
  expect_true(is.numeric(result$distance))
  expect_lt(result$distance, 1e-4)
})

test_that("shape.distance returns gamma and f2.aligned", {
  t <- seq(0, 1, length.out = 30)
  f1 <- sin(2 * pi * t)
  f2 <- cos(2 * pi * t)

  result <- shape.distance(f1, f2, argvals = t)

  expect_true(is.numeric(result$distance))
  expect_gte(result$distance, 0)
  expect_true(is.numeric(result$gamma))
  expect_equal(length(result$gamma), 30)
  expect_true(is.numeric(result$f2.aligned))
  expect_equal(length(result$f2.aligned), 30)
})

test_that("shape.distance is non-negative", {
  t <- seq(0, 1, length.out = 30)
  f1 <- sin(2 * pi * t)
  f2 <- sin(2 * pi * t^1.5)

  result <- shape.distance(f1, f2, argvals = t)
  expect_gte(result$distance, 0)
})

test_that("shape.distance with translation quotient", {
  t <- seq(0, 1, length.out = 30)
  f1 <- sin(2 * pi * t)
  # Same shape but shifted vertically
  f2 <- sin(2 * pi * t) + 10

  result <- shape.distance(f1, f2, argvals = t, quotient = "translation")

  # With translation quotient, shifted versions should have small distance
  expect_lt(result$distance, 1.0)
})

test_that("shape.distance with scale quotient", {
  t <- seq(0, 1, length.out = 30)
  f1 <- sin(2 * pi * t)
  f2 <- 5 * sin(2 * pi * t)

  result <- shape.distance(f1, f2, argvals = t, quotient = "scale")

  # With scale quotient, rescaled versions should have small distance
  expect_lt(result$distance, 1.0)
})

test_that("shape.distance default argvals", {
  f1 <- sin(2 * pi * seq(0, 1, length.out = 25))
  f2 <- cos(2 * pi * seq(0, 1, length.out = 25))

  result <- shape.distance(f1, f2)
  expect_gte(result$distance, 0)
})

# =============================================================================
# shape.mean
# =============================================================================

test_that("shape.mean returns correct structure", {
  set.seed(1)
  t <- seq(0, 1, length.out = 25)
  n <- 8
  X <- matrix(0, n, 25)
  for (i in 1:n) {
    phase <- runif(1, -0.05, 0.05)
    X[i, ] <- sin(2 * pi * (t + phase)) + rnorm(25, sd = 0.05)
  }
  fd <- fdata(X, argvals = t)

  sm <- shape.mean(fd, quotient = "reparameterization", max.iter = 5)

  expect_s3_class(sm, "shape.mean")
  expect_true(is.numeric(sm$mean))
  expect_equal(length(sm$mean), 25)
  expect_true(is.numeric(sm$mean.srsf))
  expect_true(!is.null(sm$gammas))
  expect_true(!is.null(sm$aligned.data))
  expect_true(is.numeric(sm$n.iter))
  expect_gte(sm$n.iter, 1)
  expect_true(is.logical(sm$converged))
  expect_equal(sm$quotient, "reparameterization")
})

test_that("shape.mean returns mean curve of correct length", {
  set.seed(42)
  m <- 20
  t <- seq(0, 1, length.out = m)
  n <- 6
  X <- matrix(0, n, m)
  for (i in 1:n) {
    X[i, ] <- sin(2 * pi * t) + rnorm(m, sd = 0.1)
  }
  fd <- fdata(X, argvals = t)

  sm <- shape.mean(fd, max.iter = 3)
  expect_equal(length(sm$mean), m)
})

test_that("shape.mean validates input", {
  expect_error(shape.mean("not_fdata"), "fdata")
})

test_that("shape.mean print method works", {
  set.seed(1)
  m <- 20
  t <- seq(0, 1, length.out = m)
  n <- 6
  X <- matrix(0, n, m)
  for (i in 1:n) {
    X[i, ] <- sin(2 * pi * t) + rnorm(m, sd = 0.1)
  }
  fd <- fdata(X, argvals = t)
  sm <- shape.mean(fd, max.iter = 3)
  expect_output(print(sm), "Shape Mean")
})

# =============================================================================
# shape.distance.matrix
# =============================================================================

test_that("shape.distance.matrix returns symmetric matrix", {
  set.seed(1)
  m <- 20
  t <- seq(0, 1, length.out = m)
  n <- 5
  X <- matrix(0, n, m)
  for (i in 1:n) {
    phase <- runif(1, -0.1, 0.1)
    X[i, ] <- sin(2 * pi * (t + phase)) + rnorm(m, sd = 0.1)
  }
  fd <- fdata(X, argvals = t)

  D <- shape.distance.matrix(fd)

  expect_true(is.matrix(D))
  expect_equal(nrow(D), n)
  expect_equal(ncol(D), n)
  # Symmetric
  expect_equal(D, t(D), tolerance = 1e-6)
})

test_that("shape.distance.matrix has zero diagonal", {
  set.seed(1)
  m <- 20
  t <- seq(0, 1, length.out = m)
  n <- 4
  X <- matrix(0, n, m)
  for (i in 1:n) {
    X[i, ] <- sin(2 * pi * t) + rnorm(m, sd = 0.1)
  }
  fd <- fdata(X, argvals = t)

  D <- shape.distance.matrix(fd)
  expect_true(all(diag(D) < 1e-4))
})

test_that("shape.distance.matrix entries are non-negative", {
  set.seed(42)
  m <- 20
  t <- seq(0, 1, length.out = m)
  n <- 4
  X <- matrix(0, n, m)
  for (i in 1:n) {
    X[i, ] <- sin(2 * pi * t) + rnorm(m, sd = 0.1)
  }
  fd <- fdata(X, argvals = t)

  D <- shape.distance.matrix(fd)
  expect_true(all(D >= -1e-10))
})

test_that("shape.distance.matrix with different quotients", {
  set.seed(1)
  m <- 20
  t <- seq(0, 1, length.out = m)
  n <- 4
  X <- matrix(0, n, m)
  for (i in 1:n) X[i, ] <- sin(2 * pi * t) + rnorm(m, sd = 0.1)
  fd <- fdata(X, argvals = t)

  D_reparam <- shape.distance.matrix(fd, quotient = "reparameterization")
  D_trans <- shape.distance.matrix(fd, quotient = "translation")
  D_scale <- shape.distance.matrix(fd, quotient = "scale")

  expect_true(is.matrix(D_reparam))
  expect_true(is.matrix(D_trans))
  expect_true(is.matrix(D_scale))
  expect_equal(dim(D_reparam), c(n, n))
  expect_equal(dim(D_trans), c(n, n))
  expect_equal(dim(D_scale), c(n, n))
})

test_that("shape.distance.matrix identical curves have zero distance", {
  m <- 20
  t <- seq(0, 1, length.out = m)
  X <- matrix(sin(2 * pi * t), 3, m, byrow = TRUE)
  fd <- fdata(X, argvals = t)

  D <- shape.distance.matrix(fd)
  expect_true(all(D < 1e-4))
})

test_that("shape.distance.matrix validates input", {
  expect_error(shape.distance.matrix("not_fdata"), "fdata")
})
