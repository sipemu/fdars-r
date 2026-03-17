# Tests for SPM Phase 5: Elastic SPM

# =============================================================================
# spm.elastic.phase1
# =============================================================================

test_that("spm.elastic.phase1 returns correct structure", {
  set.seed(1)
  n <- 30; m <- 25
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  for (i in 1:n) {
    phase <- runif(1, -0.05, 0.05)
    X[i, ] <- sin(2 * pi * (argvals + phase)) + rnorm(m, sd = 0.1)
  }
  fd <- fdata(X, argvals = argvals)

  chart <- spm.elastic.phase1(fd, ncomp = 3, warp.ncomp = 2, seed = 1)

  expect_s3_class(chart, "spm.elastic.chart")
  expect_true(!is.null(chart$karcher.mean))
  expect_true(is.numeric(chart$karcher.mean))
  expect_equal(length(chart$karcher.mean), m)

  # Amplitude chart fields
  expect_true(!is.null(chart$amplitude))
  expect_true(chart$amplitude$t2.ucl > 0)
  expect_true(chart$amplitude$spe.ucl > 0)
  expect_true(is.numeric(chart$amplitude$eigenvalues))
  expect_true(length(chart$amplitude$t2.phase1) > 0)
  expect_true(length(chart$amplitude$spe.phase1) > 0)

  # Phase chart enabled by default
  expect_true(chart$monitor.phase)
  expect_true(!is.null(chart$phase))
  expect_true(chart$phase$t2.ucl > 0)
  expect_true(chart$phase$spe.ucl > 0)
})

test_that("spm.elastic.phase1 without phase monitoring", {
  set.seed(1)
  n <- 25; m <- 20
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fd <- fdata(X, argvals = argvals)

  chart <- spm.elastic.phase1(fd, ncomp = 3, monitor.phase = FALSE, seed = 1)

  expect_s3_class(chart, "spm.elastic.chart")
  expect_false(chart$monitor.phase)
  expect_null(chart$phase)
  expect_true(!is.null(chart$amplitude))
})

test_that("spm.elastic.phase1 validates input", {
  expect_error(spm.elastic.phase1("not_fdata"), "fdata")
})

test_that("spm.elastic.phase1 print method works", {
  set.seed(1)
  n <- 25; m <- 20
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fd <- fdata(X, argvals = argvals)
  chart <- spm.elastic.phase1(fd, ncomp = 3, warp.ncomp = 2, seed = 1)
  expect_output(print(chart), "Elastic SPM Control Chart")
})

# =============================================================================
# spm.elastic.monitor
# =============================================================================

test_that("spm.elastic.monitor returns amplitude and phase results", {
  set.seed(1)
  n <- 30; m <- 25
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  for (i in 1:n) {
    phase <- runif(1, -0.05, 0.05)
    X[i, ] <- sin(2 * pi * (argvals + phase)) + rnorm(m, sd = 0.1)
  }
  fd <- fdata(X, argvals = argvals)
  chart <- spm.elastic.phase1(fd, ncomp = 3, warp.ncomp = 2, seed = 1)

  # New data
  set.seed(99)
  n_new <- 10
  X_new <- matrix(0, n_new, m)
  for (i in 1:n_new) {
    phase <- runif(1, -0.05, 0.05)
    X_new[i, ] <- sin(2 * pi * (argvals + phase)) + rnorm(m, sd = 0.1)
  }
  fd_new <- fdata(X_new, argvals = argvals)

  result <- spm.elastic.monitor(chart, fd_new)

  expect_true(is.list(result))
  # Amplitude results
  expect_true(!is.null(result$amplitude))
  expect_true(is.numeric(result$amplitude$t2))
  expect_equal(length(result$amplitude$t2), n_new)
  expect_true(is.logical(result$amplitude$t2.alarm))
  expect_true(is.numeric(result$amplitude$spe))
  expect_true(is.logical(result$amplitude$spe.alarm))

  # Phase results (since monitor.phase = TRUE)
  expect_true(!is.null(result$phase))
  expect_true(is.numeric(result$phase$t2))
  expect_equal(length(result$phase$t2), n_new)

  # Aligned data and warps
  expect_true(!is.null(result$aligned.data))
  expect_true(!is.null(result$warping.functions))
})

test_that("spm.elastic.monitor detects amplitude shift", {
  set.seed(1)
  n <- 30; m <- 25
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fd <- fdata(X, argvals = argvals)
  chart <- spm.elastic.phase1(fd, ncomp = 3, warp.ncomp = 2, seed = 1)

  # Shifted amplitude
  set.seed(99)
  X_shift <- matrix(rnorm(10 * m, mean = 5), 10, m)
  fd_shift <- fdata(X_shift, argvals = argvals)

  result <- spm.elastic.monitor(chart, fd_shift)
  expect_true(any(result$amplitude$t2.alarm) ||
              any(result$amplitude$spe.alarm))
})

test_that("spm.elastic.monitor without phase monitoring", {
  set.seed(1)
  n <- 25; m <- 20
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fd <- fdata(X, argvals = argvals)
  chart <- spm.elastic.phase1(fd, ncomp = 3, monitor.phase = FALSE, seed = 1)

  fd_new <- fdata(matrix(rnorm(10 * m), 10, m), argvals = argvals)
  result <- spm.elastic.monitor(chart, fd_new)

  expect_true(!is.null(result$amplitude))
  # Phase should be NULL when not monitored
  expect_null(result$phase)
})

test_that("spm.elastic.monitor validates input classes", {
  set.seed(1)
  n <- 25; m <- 20
  argvals <- seq(0, 1, length.out = m)
  fd <- fdata(matrix(rnorm(n * m), n, m), argvals = argvals)
  chart <- spm.elastic.phase1(fd, ncomp = 3, seed = 1)

  expect_error(spm.elastic.monitor("not_chart", fd), "spm.elastic.chart")
  expect_error(spm.elastic.monitor(chart, "not_fdata"), "fdata")
})
