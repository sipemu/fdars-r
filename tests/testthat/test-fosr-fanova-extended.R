# Extended tests for FoSR, fosr.fpc, fanova, predict.fosr

# =============================================================================
# Helper
# =============================================================================

.test_fosr_data <- function(n = 50, m = 20, p = 2, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(0, n, m)
  for (i in 1:n) {
    Y[i, ] <- 3 * sin(2 * pi * t) + X[i, 1] * cos(2 * pi * t) +
              X[i, 2] * t + rnorm(m, sd = 0.3)
  }
  fd <- fdata(Y, argvals = t)
  list(fd = fd, X = X)
}

# =============================================================================
# fosr with multiple predictors
# =============================================================================

test_that("fosr returns correct structure with multiple predictors", {
  d <- .test_fosr_data()
  fit <- fosr(d$fd, d$X)

  expect_s3_class(fit, "fosr")
  expect_s3_class(fit$intercept, "fdata")
  expect_s3_class(fit$beta, "fdata")
  expect_s3_class(fit$fitted, "fdata")
  expect_s3_class(fit$residuals, "fdata")
  expect_true(fit$r.squared >= 0 && fit$r.squared <= 1)
  expect_equal(nrow(fit$fitted$data), 50)
  expect_equal(nrow(fit$beta$data), 2)
})

test_that("fosr with lambda penalization", {
  d <- .test_fosr_data()
  fit_nopenalty <- fosr(d$fd, d$X, lambda = 0)
  fit_penalty <- fosr(d$fd, d$X, lambda = 10)

  expect_s3_class(fit_penalty, "fosr")
  # Penalized version should have smoother coefficients (smaller norm)
  beta_norm_0 <- sum(fit_nopenalty$beta$data^2)
  beta_norm_10 <- sum(fit_penalty$beta$data^2)
  expect_lt(beta_norm_10, beta_norm_0)
})

test_that("fosr residuals + fitted = data", {
  d <- .test_fosr_data()
  fit <- fosr(d$fd, d$X)

  reconstructed <- fit$fitted$data + fit$residuals$data
  expect_equal(reconstructed, d$fd$data, tolerance = 1e-8)
})

test_that("fosr r.squared.t is valid", {
  d <- .test_fosr_data()
  fit <- fosr(d$fd, d$X)

  expect_equal(length(fit$r.squared.t), ncol(d$fd$data))
  expect_true(all(is.finite(fit$r.squared.t)))
})

# =============================================================================
# predict.fosr
# =============================================================================

test_that("predict.fosr with NULL returns fitted", {
  d <- .test_fosr_data()
  fit <- fosr(d$fd, d$X)

  pred <- predict(fit)
  expect_s3_class(pred, "fdata")
  expect_equal(pred$data, fit$fitted$data)
})

test_that("predict.fosr with new data returns correct dims", {
  d <- .test_fosr_data()
  fit <- fosr(d$fd, d$X)

  X_new <- matrix(rnorm(10), 5, 2)
  pred <- predict(fit, X_new)
  expect_s3_class(pred, "fdata")
  expect_equal(nrow(pred$data), 5)
  expect_equal(ncol(pred$data), ncol(d$fd$data))
})

# =============================================================================
# fosr.fpc
# =============================================================================

test_that("fosr.fpc returns correct structure", {
  d <- .test_fosr_data()
  fit <- fosr.fpc(d$fd, d$X, ncomp = 3)

  expect_s3_class(fit, "fosr")
  expect_true(fit$r.squared >= 0 && fit$r.squared <= 1)
  expect_equal(fit$ncomp, 3)
})

test_that("fosr.fpc with varying ncomp", {
  d <- .test_fosr_data()
  fit2 <- fosr.fpc(d$fd, d$X, ncomp = 2)
  fit5 <- fosr.fpc(d$fd, d$X, ncomp = 5)

  expect_equal(fit2$ncomp, 2)
  expect_equal(fit5$ncomp, 5)
  # More components should generally give higher R²
  expect_gte(fit5$r.squared, fit2$r.squared - 0.05)
})

# =============================================================================
# fanova
# =============================================================================

test_that("fanova detects group differences", {
  set.seed(42)
  m <- 20
  t <- seq(0, 1, length.out = m)
  n <- 30
  X <- matrix(0, n, m)
  groups <- rep(1:3, each = 10)
  for (i in 1:n) {
    g <- groups[i]
    X[i, ] <- g * sin(2 * pi * t) + rnorm(m, sd = 0.3)
  }
  fd <- fdata(X, argvals = t)

  result <- fanova(fd, groups)

  expect_s3_class(result, "fanova")
  expect_true(result$p.value < 0.05)  # Strong signal should be significant
  expect_equal(result$n.groups, 3)
  expect_s3_class(result$group.means, "fdata")
})

test_that("fanova p-value is in [0, 1]", {
  set.seed(42)
  m <- 20
  X <- matrix(rnorm(40 * m), 40, m)
  fd <- fdata(X, argvals = seq(0, 1, length.out = m))
  groups <- rep(1:2, each = 20)

  result <- fanova(fd, groups, n.perm = 50)
  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("fanova validates mismatched groups", {
  fd <- fdata(matrix(rnorm(100), 10, 10))
  expect_error(fanova(fd, 1:5), "groups")
})

# =============================================================================
# fosr.2d extended
# =============================================================================

test_that("fosr.2d with lambda penalization", {
  set.seed(42)
  n <- 30
  m_s <- 8; m_t <- 6
  m <- m_s * m_t
  Y <- matrix(rnorm(n * m), n, m)
  fd <- fdata(Y, argvals = seq_len(m))
  X <- cbind(rnorm(n))

  fit <- fosr.2d(fd, X, seq(0, 1, length.out = m_s),
                 seq(0, 1, length.out = m_t), lambda.s = 1.0, lambda.t = 1.0)

  expect_s3_class(fit, "fosr.2d")
  expect_true(is.numeric(fit$r.squared))
})

test_that("fosr.2d with multiple predictors", {
  set.seed(42)
  n <- 30
  m_s <- 8; m_t <- 6
  m <- m_s * m_t
  Y <- matrix(rnorm(n * m), n, m)
  fd <- fdata(Y, argvals = seq_len(m))
  X <- cbind(rnorm(n), rnorm(n))

  fit <- fosr.2d(fd, X, seq(0, 1, length.out = m_s),
                 seq(0, 1, length.out = m_t))

  expect_s3_class(fit, "fosr.2d")
  expect_equal(nrow(fit$beta$data), 2)
})

# =============================================================================
# plot methods
# =============================================================================

test_that("plot.fosr works", {
  d <- .test_fosr_data()
  fit <- fosr(d$fd, d$X)

  pdf(NULL)
  expect_no_error(plot(fit))
  dev.off()
})

test_that("plot.fanova works", {
  set.seed(42)
  m <- 20
  X <- matrix(rnorm(30 * m), 30, m)
  fd <- fdata(X, argvals = seq(0, 1, length.out = m))
  groups <- rep(1:3, each = 10)
  result <- fanova(fd, groups)

  pdf(NULL)
  expect_no_error(plot(result))
  dev.off()
})
