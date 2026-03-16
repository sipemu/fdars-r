# Tests for alignment utility functions

# =============================================================================
# elastic.pair
# =============================================================================

test_that("elastic.pair returns correct structure", {
  t <- seq(0, 1, length.out = 50)
  f1 <- fdata(matrix(sin(2 * pi * t), 1, 50), argvals = t)
  f2 <- fdata(matrix(sin(2 * pi * t^1.5), 1, 50), argvals = t)

  result <- elastic.pair(f1, f2)

  expect_true(is.list(result))
  expect_s3_class(result$f.aligned, "fdata")
  expect_true(is.numeric(result$gamma))
  expect_equal(length(result$gamma), 50)
  expect_true(is.numeric(result$distance))
  expect_true(result$distance >= 0)
})

test_that("elastic.pair aligns identical curves with zero distance", {
  t <- seq(0, 1, length.out = 50)
  f1 <- fdata(matrix(sin(2 * pi * t), 1, 50), argvals = t)

  result <- elastic.pair(f1, f1)

  expect_true(result$distance < 1e-6)
})

test_that("elastic.pair reduces distance after alignment", {
  t <- seq(0, 1, length.out = 50)
  f1 <- fdata(matrix(sin(2 * pi * t), 1, 50), argvals = t)
  f2 <- fdata(matrix(sin(2 * pi * t^1.5), 1, 50), argvals = t)

  result <- elastic.pair(f1, f2)

  # Pre-alignment L2 distance should be larger
  pre_dist <- sqrt(sum((f1$data[1, ] - f2$data[1, ])^2) / 50)
  expect_true(result$distance <= pre_dist + 1e-6)
})

test_that("elastic.pair validates input", {
  expect_error(elastic.pair("not_fdata", "not_fdata"), "fdata")
})

# =============================================================================
# warp.compose
# =============================================================================

test_that("warp.compose returns correct length", {
  t <- seq(0, 1, length.out = 20)
  gamma1 <- t^1.5
  gamma2 <- t^0.8
  fd <- fdata(matrix(rnorm(20), 1, 20), argvals = t)

  result <- warp.compose(gamma1, gamma2, fd)

  expect_equal(length(result), 20)
  expect_true(is.numeric(result))
})

test_that("warp.compose with identity preserves warp", {
  t <- seq(0, 1, length.out = 20)
  gamma <- t^1.5
  identity <- t

  result <- warp.compose(gamma, identity, t)

  expect_equal(result, gamma, tolerance = 1e-6)
})

test_that("warp.compose accepts argvals vector", {
  t <- seq(0, 1, length.out = 20)
  gamma1 <- t^1.3
  gamma2 <- t^0.9

  result <- warp.compose(gamma1, gamma2, t)
  expect_equal(length(result), 20)
})

# =============================================================================
# warp.complexity
# =============================================================================

test_that("warp.complexity of identity is zero", {
  t <- seq(0, 1, length.out = 30)
  expect_equal(warp.complexity(t, t), 0, tolerance = 1e-10)
})

test_that("warp.complexity is positive for non-identity", {
  t <- seq(0, 1, length.out = 30)
  gamma <- t^2
  expect_true(warp.complexity(gamma, t) > 0)
})

test_that("warp.complexity increases with distortion", {
  t <- seq(0, 1, length.out = 30)
  c1 <- warp.complexity(t^1.5, t)
  c2 <- warp.complexity(t^3, t)
  expect_true(c2 > c1)
})

# =============================================================================
# warp.smoothness
# =============================================================================

test_that("warp.smoothness returns non-negative value", {
  t <- seq(0, 1, length.out = 30)
  expect_true(warp.smoothness(t^2, t) >= 0)
})

test_that("warp.smoothness of identity is near zero", {
  t <- seq(0, 1, length.out = 30)
  expect_true(warp.smoothness(t, t) < 1e-6)
})
