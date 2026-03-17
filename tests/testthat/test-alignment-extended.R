# Extended alignment tests: pairwise consistency, metric.elastic, srsf.reparameterize

# =============================================================================
# Helper
# =============================================================================

.test_align_data <- function(n = 10, m = 50, seed = 42) {
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
# alignment.pairwise.consistency
# =============================================================================

test_that("alignment.pairwise.consistency returns scalar in [0, 1]", {
  fd <- .test_align_data(n = 6, m = 30)
  result <- alignment.pairwise.consistency(fd, lambda = 0)

  expect_true(is.numeric(result))
  expect_equal(length(result), 1)
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("alignment.pairwise.consistency returns finite value for similar curves", {
  set.seed(42)
  t <- seq(0, 1, length.out = 30)
  X <- matrix(0, 5, 30)
  for (i in 1:5) {
    X[i, ] <- sin(2 * pi * t) + rnorm(30, sd = 0.01)
  }
  fd <- fdata(X, argvals = t)

  result <- alignment.pairwise.consistency(fd)
  expect_true(is.finite(result))
  expect_gte(result, 0)
  expect_lte(result, 1)
})

test_that("alignment.pairwise.consistency with max.triplets", {
  fd <- .test_align_data(n = 8, m = 30)

  result_all <- alignment.pairwise.consistency(fd, max.triplets = 0)
  result_sub <- alignment.pairwise.consistency(fd, max.triplets = 10)

  expect_true(is.numeric(result_sub))
  expect_gte(result_sub, 0)
  expect_lte(result_sub, 1)
})

test_that("alignment.pairwise.consistency validates input", {
  expect_error(alignment.pairwise.consistency("not_fdata"), "fdata")
})

# =============================================================================
# metric.elastic
# =============================================================================

test_that("metric.elastic returns distance matrix", {
  fd <- .test_align_data(n = 5, m = 30)
  result <- metric.elastic(fd)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 5)
  # Distance matrix should be symmetric
  expect_equal(result, t(result), tolerance = 1e-8)
  # Diagonal should be zero (or near-zero)
  expect_true(all(diag(result) < 1e-6))
  # Off-diagonal should be non-negative
  expect_true(all(result >= -1e-10))
})

test_that("metric.elastic with reference data", {
  fd <- .test_align_data(n = 5, m = 30)
  fd_ref <- .test_align_data(n = 3, m = 30, seed = 99)

  result <- metric.elastic(fd, fd_ref)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 5)
  expect_equal(ncol(result), 3)
})

test_that("metric.elastic identical curves have zero distance", {
  t <- seq(0, 1, length.out = 30)
  X <- matrix(sin(2 * pi * t), 3, 30, byrow = TRUE)
  fd <- fdata(X, argvals = t)

  result <- metric.elastic(fd)
  expect_true(all(result < 1e-6))
})

# =============================================================================
# srsf.reparameterize
# =============================================================================

test_that("srsf.reparameterize returns fdata", {
  t <- seq(0, 1, length.out = 30)
  fd <- fdata(matrix(sin(2 * pi * t), 1, 30), argvals = t)
  gamma <- t^2  # simple warping function

  result <- srsf.reparameterize(fd, gamma)

  expect_s3_class(result, "fdata")
  expect_equal(nrow(result$data), 1)
  expect_equal(ncol(result$data), 30)
  expect_equal(result$argvals, t)
})

test_that("srsf.reparameterize with identity warp preserves curve", {
  t <- seq(0, 1, length.out = 30)
  fd <- fdata(matrix(sin(2 * pi * t), 1, 30), argvals = t)
  gamma <- t  # identity warping

  result <- srsf.reparameterize(fd, gamma)
  expect_equal(as.vector(result$data), sin(2 * pi * t), tolerance = 1e-4)
})

test_that("srsf.reparameterize validates input", {
  expect_error(srsf.reparameterize("not_fdata", 1:10), "fdata")
})
