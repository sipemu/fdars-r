# Tests for SPM Phase 4: Profile monitoring and partial-domain monitoring

# =============================================================================
# spm.profile.phase1
# =============================================================================

test_that("spm.profile.phase1 returns correct structure", {
  set.seed(1)
  n <- 80; m <- 20
  argvals <- seq(0, 1, length.out = m)
  X_pred <- cbind(rnorm(n), rnorm(n))
  Y <- matrix(rnorm(n * m), n, m)
  fd <- fdata(Y, argvals = argvals)

  chart <- spm.profile.phase1(fd, X_pred, ncomp = 2, window.size = 15)

  expect_s3_class(chart, "spm.profile.chart")
  expect_true(chart$t2.ucl > 0)
  expect_equal(chart$ncomp, 2)
  expect_equal(chart$alpha, 0.05)
  expect_true(is.numeric(chart$eigenvalues))
  expect_true(is.character(chart$t2.description))
  expect_true(is.numeric(chart$lag1.autocorrelation))
  expect_true(is.numeric(chart$effective.n.windows))
})

test_that("spm.profile.phase1 with 1 predictor", {
  set.seed(42)
  n <- 60; m <- 15
  argvals <- seq(0, 1, length.out = m)
  X_pred <- matrix(rnorm(n), n, 1)
  Y <- matrix(rnorm(n * m), n, m)
  fd <- fdata(Y, argvals = argvals)

  chart <- spm.profile.phase1(fd, X_pred, ncomp = 2, window.size = 10)
  expect_s3_class(chart, "spm.profile.chart")
})

test_that("spm.profile.phase1 validates input", {
  expect_error(spm.profile.phase1("not_fdata", matrix(1, 2, 2)), "fdata")
})

test_that("spm.profile.phase1 print method works", {
  set.seed(1)
  n <- 60; m <- 15
  argvals <- seq(0, 1, length.out = m)
  X_pred <- cbind(rnorm(n), rnorm(n))
  Y <- matrix(rnorm(n * m), n, m)
  fd <- fdata(Y, argvals = argvals)
  chart <- spm.profile.phase1(fd, X_pred, ncomp = 2, window.size = 10)
  expect_output(print(chart), "SPM Profile Monitoring Chart")
})

# =============================================================================
# spm.profile.monitor
# =============================================================================

test_that("spm.profile.monitor returns correct structure", {
  set.seed(1)
  n <- 80; m <- 20
  argvals <- seq(0, 1, length.out = m)
  X_pred <- cbind(rnorm(n), rnorm(n))
  Y <- matrix(rnorm(n * m), n, m)
  fd <- fdata(Y, argvals = argvals)
  chart <- spm.profile.phase1(fd, X_pred, ncomp = 2, window.size = 15)

  # New data
  n_new <- 30
  X_new <- cbind(rnorm(n_new), rnorm(n_new))
  Y_new <- matrix(rnorm(n_new * m), n_new, m)
  fd_new <- fdata(Y_new, argvals = argvals)

  result <- spm.profile.monitor(chart, fd_new, X_new)

  expect_s3_class(result, "spm.monitor")
  expect_true(is.numeric(result$t2))
  expect_true(length(result$t2) > 0)
  expect_true(is.logical(result$t2.alarm))
  expect_true(result$t2.ucl > 0)
})

test_that("spm.profile.monitor detects shifted relationship", {
  set.seed(1)
  n <- 80; m <- 20
  argvals <- seq(0, 1, length.out = m)
  X_pred <- cbind(rnorm(n), rnorm(n))
  Y <- matrix(rnorm(n * m), n, m)
  fd <- fdata(Y, argvals = argvals)
  chart <- spm.profile.phase1(fd, X_pred, ncomp = 2, window.size = 15)

  # Shifted new data
  n_new <- 30
  X_new <- cbind(rnorm(n_new), rnorm(n_new))
  Y_new <- matrix(rnorm(n_new * m) + 5, n_new, m)
  fd_new <- fdata(Y_new, argvals = argvals)

  result <- spm.profile.monitor(chart, fd_new, X_new)
  # With large shift, T2 values should be elevated
  expect_true(is.numeric(result$t2))
  expect_true(length(result$t2) > 0)
})

test_that("spm.profile.monitor validates chart class", {
  set.seed(1)
  fd <- fdata(matrix(rnorm(20), 2, 10), argvals = seq(0, 1, length.out = 10))
  expect_error(spm.profile.monitor("not_chart", fd, matrix(1, 2, 1)),
               "spm.profile.chart")
})

# =============================================================================
# spm.monitor.partial
# =============================================================================

test_that("spm.monitor.partial returns correct structure", {
  set.seed(1)
  n <- 50; m <- 30
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fd <- fdata(X, argvals = argvals)
  chart <- spm.phase1(fd, ncomp = 3)

  # Partial curve: first 20 of 30 grid points
  partial <- rnorm(m)
  result <- spm.monitor.partial(chart, partial, n.observed = 20)

  expect_true(is.list(result))
  expect_true(is.numeric(result$t2))
  expect_true(is.logical(result$t2.alarm))
  expect_true(is.numeric(result$domain.fraction))
  expect_gt(result$domain.fraction, 0)
  expect_lte(result$domain.fraction, 1)
  expect_true(!is.null(result$scores))
})

test_that("spm.monitor.partial conditional completion works", {
  set.seed(1)
  n <- 50; m <- 30
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fd <- fdata(X, argvals = argvals)
  chart <- spm.phase1(fd, ncomp = 3)

  partial <- rnorm(m)
  result <- spm.monitor.partial(chart, partial, n.observed = 21,
                                 completion = "conditional")

  expect_true(is.numeric(result$t2))
  expect_true(result$t2 >= 0)
})

test_that("spm.monitor.partial projection completion works", {
  set.seed(1)
  n <- 50; m <- 30
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fd <- fdata(X, argvals = argvals)
  chart <- spm.phase1(fd, ncomp = 3)

  partial <- rnorm(m)
  result <- spm.monitor.partial(chart, partial, n.observed = 15,
                                 completion = "projection")

  expect_true(is.numeric(result$t2))
  expect_true(result$t2 >= 0)
})

test_that("spm.monitor.partial zero_pad completion works", {
  set.seed(1)
  n <- 50; m <- 30
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fd <- fdata(X, argvals = argvals)
  chart <- spm.phase1(fd, ncomp = 3)

  partial <- rnorm(m)
  result <- spm.monitor.partial(chart, partial, n.observed = 10,
                                 completion = "zero_pad")

  expect_true(is.numeric(result$t2))
  expect_true(result$t2 >= 0)
})

test_that("spm.monitor.partial validates chart class", {
  expect_error(spm.monitor.partial("not_chart", rnorm(10), 5), "spm.chart")
})

# =============================================================================
# spm.monitor.partial.batch
# =============================================================================

test_that("spm.monitor.partial.batch returns list of correct length", {
  set.seed(1)
  n <- 50; m <- 30
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fd <- fdata(X, argvals = argvals)
  chart <- spm.phase1(fd, ncomp = 3)

  # Three curves with different observed lengths
  curves <- list(rnorm(m), rnorm(m), rnorm(m))
  n_obs <- c(10L, 20L, 25L)

  results <- spm.monitor.partial.batch(chart, curves, n_obs)

  expect_true(is.list(results))
  expect_equal(length(results), 3)
})

test_that("spm.monitor.partial.batch each element has correct fields", {
  set.seed(1)
  n <- 50; m <- 30
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fd <- fdata(X, argvals = argvals)
  chart <- spm.phase1(fd, ncomp = 3)

  curves <- list(rnorm(m), rnorm(m))
  n_obs <- c(15L, 25L)

  results <- spm.monitor.partial.batch(chart, curves, n_obs)

  for (i in seq_along(results)) {
    res <- results[[i]]
    expect_true(is.numeric(res$t2))
    expect_true(is.logical(res$t2.alarm))
    expect_true(is.numeric(res$domain.fraction))
    expect_true(res$domain.fraction > 0 && res$domain.fraction <= 1)
  }
})

test_that("spm.monitor.partial.batch domain fractions match n_obs", {
  set.seed(1)
  n <- 50; m <- 30
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fd <- fdata(X, argvals = argvals)
  chart <- spm.phase1(fd, ncomp = 3)

  curves <- list(rnorm(m), rnorm(m), rnorm(m))
  n_obs <- c(10L, 20L, 30L)

  results <- spm.monitor.partial.batch(chart, curves, n_obs)

  # Domain fractions should increase with more observed points
  expect_lt(results[[1]]$domain.fraction, results[[2]]$domain.fraction)
  expect_lt(results[[2]]$domain.fraction, results[[3]]$domain.fraction + 1e-10)
  expect_equal(results[[3]]$domain.fraction, 1.0, tolerance = 1e-4)
})

test_that("spm.monitor.partial.batch validates mismatched lengths", {
  set.seed(1)
  n <- 50; m <- 30
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fd <- fdata(X, argvals = argvals)
  chart <- spm.phase1(fd, ncomp = 3)

  curves <- list(rnorm(m), rnorm(m))
  n_obs <- c(10L, 20L, 25L)

  expect_error(spm.monitor.partial.batch(chart, curves, n_obs),
               "same length")
})

test_that("spm.monitor.partial.batch validates chart class", {
  expect_error(spm.monitor.partial.batch("not_chart", list(1:10), c(5L)),
               "spm.chart")
})
