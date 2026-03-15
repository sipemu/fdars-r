# Extended tests for tolerance band methods

# Helper: create smooth functional data
.test_tol_data <- function(n = 50, m = 20, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  for (i in 1:n) {
    X[i, ] <- sin(2 * pi * t) + rnorm(1, 0, 0.3) * cos(4 * pi * t) +
              rnorm(m, sd = 0.1)
  }
  fdata(X, argvals = t)
}

# =============================================================================
# Method coverage
# =============================================================================

test_that("tolerance.band exponential method works", {
  fd <- .test_tol_data()
  band <- tolerance.band(fd, method = "exponential", nb = 50, seed = 42)
  expect_s3_class(band, "tolerance.band")
  expect_equal(band$method, "exponential")
  expect_equal(length(band$lower), 20)
  expect_equal(length(band$upper), 20)
  expect_true(all(band$upper >= band$lower))
})

test_that("tolerance.band elastic method works", {
  fd <- .test_tol_data()
  band <- tolerance.band(fd, method = "elastic", nb = 50, ncomp = 2,
                          max.iter = 3, seed = 42)
  expect_s3_class(band, "tolerance.band")
  expect_equal(band$method, "elastic")
  expect_true(all(band$upper >= band$lower))
})

test_that("tolerance.band phase method works", {
  fd <- .test_tol_data()
  band <- tolerance.band(fd, method = "phase", nb = 50, ncomp = 2,
                          max.iter = 3, seed = 42)
  expect_s3_class(band, "tolerance.band")
  expect_equal(band$method, "phase")
  expect_true(!is.null(band$gamma.lower) || !is.null(band$phase.lower))
})

test_that("tolerance.band elastic.config method works", {
  fd <- .test_tol_data()
  band <- tolerance.band(fd, method = "elastic.config", nb = 50, ncomp = 2,
                          ncomp.phase = 2, max.iter = 3, seed = 42)
  expect_s3_class(band, "tolerance.band")
  expect_equal(band$method, "elastic.config")
  expect_true(all(band$upper >= band$lower))
})

# =============================================================================
# Coverage parameter
# =============================================================================

test_that("tolerance.band respects coverage parameter", {
  fd <- .test_tol_data()
  band_90 <- tolerance.band(fd, method = "fpca", coverage = 0.90, nb = 50, seed = 42)
  band_99 <- tolerance.band(fd, method = "fpca", coverage = 0.99, nb = 50, seed = 42)

  expect_equal(band_90$coverage, 0.90)
  expect_equal(band_99$coverage, 0.99)
  # Higher coverage should give wider bands
  width_90 <- mean(band_90$upper - band_90$lower)
  width_99 <- mean(band_99$upper - band_99$lower)
  expect_true(width_99 >= width_90)
})

# =============================================================================
# Plot method
# =============================================================================

test_that("plot.tolerance.band runs without error", {
  fd <- .test_tol_data()
  band <- tolerance.band(fd, method = "fpca", nb = 50, seed = 42)
  pdf(NULL)
  expect_no_error(plot(band))
  dev.off()
})

# =============================================================================
# Edge cases
# =============================================================================

test_that("tolerance.band with few observations works", {
  fd <- .test_tol_data(n = 10, m = 10)
  band <- tolerance.band(fd, method = "fpca", ncomp = 2, nb = 20, seed = 42)
  expect_s3_class(band, "tolerance.band")
})
