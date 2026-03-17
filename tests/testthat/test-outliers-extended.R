# Extended outlier tests: outliers.lrt.dist, edge cases

# =============================================================================
# outliers.lrt.dist
# =============================================================================

test_that("outliers.lrt.dist returns threshold and distribution", {
  set.seed(42)
  m <- 20
  t <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(40 * m), 40, m)
  fd <- fdata(X, argvals = t)

  result <- outliers.lrt.dist(fd, nb = 50, seed = 42)

  expect_true(is.list(result))
  expect_true(is.numeric(result$threshold))
  expect_true(result$threshold > 0)
  expect_true(is.numeric(result$boot.distribution))
  expect_equal(length(result$boot.distribution), 50)
  expect_true(all(result$boot.distribution >= 0))
})

test_that("outliers.lrt.dist threshold increases with percentile", {
  set.seed(42)
  m <- 20
  fd <- fdata(matrix(rnorm(40 * m), 40, m), argvals = seq(0, 1, length.out = m))

  r_90 <- outliers.lrt.dist(fd, nb = 50, percentile = 0.90, seed = 42)
  r_99 <- outliers.lrt.dist(fd, nb = 50, percentile = 0.99, seed = 42)

  expect_lte(r_90$threshold, r_99$threshold)
})

test_that("outliers.lrt.dist validates input", {
  expect_error(outliers.lrt.dist("not_fdata"), "fdata")
})

test_that("outliers.lrt.dist rejects invalid percentile", {
  fd <- fdata(matrix(rnorm(100), 10, 10))
  expect_error(outliers.lrt.dist(fd, percentile = 0), "percentile")
  expect_error(outliers.lrt.dist(fd, percentile = 1), "percentile")
})

test_that("outliers.lrt.dist with custom trim and smo", {
  set.seed(42)
  m <- 20
  fd <- fdata(matrix(rnorm(40 * m), 40, m), argvals = seq(0, 1, length.out = m))

  result <- outliers.lrt.dist(fd, nb = 30, smo = 0.1, trim = 0.2, seed = 42)
  expect_true(result$threshold > 0)
  expect_equal(length(result$boot.distribution), 30)
})

# =============================================================================
# outliers.lrt extended
# =============================================================================

test_that("outliers.lrt detects extreme outliers", {
  set.seed(42)
  m <- 20
  t <- seq(0, 1, length.out = m)
  n <- 40
  X <- matrix(rnorm(n * m), n, m)
  # Add one extreme outlier
  X[1, ] <- 10 + rnorm(m, sd = 0.1)
  fd <- fdata(X, argvals = t)

  result <- outliers.lrt(fd, nb = 50, seed = 42)

  expect_s3_class(result, "outliers.fdata")
  # Curve 1 should be flagged as outlier
  expect_true(1 %in% result$outliers)
})

test_that("outliers.thres.lrt returns threshold", {
  set.seed(42)
  m <- 20
  fd <- fdata(matrix(rnorm(40 * m), 40, m), argvals = seq(0, 1, length.out = m))

  thresh <- outliers.thres.lrt(fd, nb = 50, seed = 42)
  expect_true(is.numeric(thresh))
  expect_true(thresh > 0)
})

# =============================================================================
# outliers.depth.pond extended
# =============================================================================

test_that("outliers.depth.pond returns correct structure", {
  set.seed(42)
  m <- 20
  n <- 30
  X <- matrix(rnorm(n * m), n, m)
  # Add outliers
  X[1, ] <- 5 + rnorm(m, sd = 0.1)
  fd <- fdata(X, argvals = seq(0, 1, length.out = m))

  result <- outliers.depth.pond(fd)

  expect_s3_class(result, "outliers.fdata")
  expect_true(is.integer(result$outliers) || is.numeric(result$outliers))
})

test_that("outliers.depth.trim returns correct structure", {
  set.seed(42)
  m <- 20
  n <- 30
  X <- matrix(rnorm(n * m), n, m)
  X[1, ] <- 5 + rnorm(m, sd = 0.1)
  fd <- fdata(X, argvals = seq(0, 1, length.out = m))

  result <- outliers.depth.trim(fd)

  expect_s3_class(result, "outliers.fdata")
  expect_true(is.integer(result$outliers) || is.numeric(result$outliers))
})
