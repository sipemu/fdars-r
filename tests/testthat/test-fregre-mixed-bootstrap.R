# Tests for fregre.np.mixed and fregre.bootstrap.ci

# =============================================================================
# Helper
# =============================================================================

.test_mixed_data <- function(n = 50, m = 20, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  z <- matrix(rnorm(n * 2), n, 2)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(m, sd = 0.1)
    y[i] <- 2 * amp + 0.5 * z[i, 1] - 0.3 * z[i, 2] + rnorm(1, sd = 0.2)
  }
  fd <- fdata(X, argvals = t)
  list(fd = fd, y = y, z = z)
}

# =============================================================================
# fregre.np.mixed
# =============================================================================

test_that("fregre.np.mixed returns correct structure", {
  d <- .test_mixed_data()
  fit <- fregre.np.mixed(d$fd, d$y, scalar.covariates = d$z)

  expect_s3_class(fit, "fregre.np")
  expect_equal(length(fit$fitted.values), 50)
  expect_equal(length(fit$residuals), 50)
  expect_true(is.numeric(fit$r.squared))
  expect_true(fit$h.func > 0)
  expect_true(fit$h.scalar > 0)
})

test_that("fregre.np.mixed works without scalar covariates", {
  d <- .test_mixed_data()
  fit <- fregre.np.mixed(d$fd, d$y)

  expect_s3_class(fit, "fregre.np")
  expect_equal(length(fit$fitted.values), 50)
  expect_true(is.numeric(fit$r.squared))
})

test_that("fregre.np.mixed fitted + residuals = y", {
  d <- .test_mixed_data()
  fit <- fregre.np.mixed(d$fd, d$y, scalar.covariates = d$z)

  expect_equal(fit$fitted.values + fit$residuals, d$y, tolerance = 1e-10)
})

test_that("fregre.np.mixed validates fdata input", {
  expect_error(fregre.np.mixed("not_fdata", 1:10), "fdata")
})

test_that("fregre.np.mixed validates y length", {
  d <- .test_mixed_data()
  expect_error(fregre.np.mixed(d$fd, 1:5), "Length of y")
})

test_that("fregre.np.mixed with explicit bandwidths", {
  d <- .test_mixed_data()
  fit <- fregre.np.mixed(d$fd, d$y, scalar.covariates = d$z,
                          h.func = 0.5, h.scalar = 0.5)

  expect_s3_class(fit, "fregre.np")
  expect_equal(length(fit$fitted.values), 50)
})

# =============================================================================
# fregre.bootstrap.ci with fregre.lm
# =============================================================================

test_that("fregre.bootstrap.ci returns correct structure for lm", {
  d <- .test_mixed_data()
  model <- fregre.lm(d$fd, d$y, ncomp = 3)
  ci <- fregre.bootstrap.ci(model, n.boot = 30, seed = 42)

  expect_s3_class(ci, "fregre.bootstrap.ci")
  expect_equal(length(ci$lower), ncol(d$fd$data))
  expect_equal(length(ci$upper), ncol(d$fd$data))
  expect_equal(length(ci$center), ncol(d$fd$data))
  expect_equal(length(ci$sim.lower), ncol(d$fd$data))
  expect_equal(length(ci$sim.upper), ncol(d$fd$data))
  expect_true(ci$n.boot.success > 0)
  expect_equal(ci$alpha, 0.05)
  expect_equal(ci$model.class, "fregre.lm")
})

test_that("fregre.bootstrap.ci lower <= center <= upper", {
  d <- .test_mixed_data()
  model <- fregre.lm(d$fd, d$y, ncomp = 3)
  ci <- fregre.bootstrap.ci(model, n.boot = 30, seed = 42)

  expect_true(all(ci$lower <= ci$center + 1e-10))
  expect_true(all(ci$center <= ci$upper + 1e-10))
  expect_true(all(ci$sim.lower <= ci$sim.upper + 1e-10))
})

test_that("fregre.bootstrap.ci with smaller alpha gives wider intervals", {
  d <- .test_mixed_data()
  model <- fregre.lm(d$fd, d$y, ncomp = 3)
  ci_90 <- fregre.bootstrap.ci(model, n.boot = 50, alpha = 0.10, seed = 42)
  ci_99 <- fregre.bootstrap.ci(model, n.boot = 50, alpha = 0.01, seed = 42)

  width_90 <- mean(ci_90$upper - ci_90$lower)
  width_99 <- mean(ci_99$upper - ci_99$lower)
  expect_true(width_99 >= width_90)
})

test_that("fregre.bootstrap.ci with logistic model", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 40
  X <- matrix(0, n, 20)
  for (i in 1:n) X[i, ] <- sin(2*pi*t) * (2*(i > 20) - 1) + rnorm(20, sd = 0.2)
  fd <- fdata(X, argvals = t)
  y <- as.numeric(1:n > 20)

  model <- functional.logistic(fd, y, ncomp = 3)
  ci <- fregre.bootstrap.ci(model, n.boot = 20, seed = 42)

  expect_s3_class(ci, "fregre.bootstrap.ci")
  expect_equal(ci$model.class, "fregre.logistic")
  expect_true(ci$n.boot.success > 0)
})

test_that("fregre.bootstrap.ci rejects invalid model class", {
  expect_error(fregre.bootstrap.ci(list(a = 1)),
               "fregre.lm.*fregre.logistic")
})

test_that("print.fregre.bootstrap.ci works", {
  d <- .test_mixed_data()
  model <- fregre.lm(d$fd, d$y, ncomp = 3)
  ci <- fregre.bootstrap.ci(model, n.boot = 20, seed = 42)

  expect_output(print(ci), "Bootstrap Confidence Intervals")
  expect_output(print(ci), "Alpha")
})

test_that("plot.fregre.bootstrap.ci works", {
  d <- .test_mixed_data()
  model <- fregre.lm(d$fd, d$y, ncomp = 3)
  ci <- fregre.bootstrap.ci(model, n.boot = 20, seed = 42)

  pdf(NULL)
  expect_no_error(plot(ci))
  expect_no_error(plot(ci, simultaneous = TRUE))
  dev.off()
})
