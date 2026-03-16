# Tests for Statistical Process Monitoring (SPM) functions

# Helper: create in-control functional data
.test_spm_data <- function(n = 50, m = 30, seed = 42) {
  set.seed(seed)
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fdata(X, argvals = argvals)
}

# =============================================================================
# spm.phase1
# =============================================================================

test_that("spm.phase1 returns correct structure", {
  fd <- .test_spm_data()
  chart <- spm.phase1(fd, ncomp = 3)

  expect_s3_class(chart, "spm.chart")
  expect_equal(chart$ncomp, 3L)
  expect_equal(chart$alpha, 0.05)
  expect_true(chart$t2.ucl > 0)
  expect_true(chart$spe.ucl > 0)
  expect_true(length(chart$t2.phase1) > 0)
  expect_true(length(chart$spe.phase1) > 0)
  expect_equal(length(chart$t2.phase1), length(chart$spe.phase1))
  expect_equal(length(chart$eigenvalues), 3)
  expect_true(all(chart$eigenvalues > 0))
  expect_true(is.character(chart$t2.description))
  expect_true(is.character(chart$spe.description))
})

test_that("spm.phase1 respects ncomp parameter", {
  fd <- .test_spm_data()
  chart2 <- spm.phase1(fd, ncomp = 2)
  chart5 <- spm.phase1(fd, ncomp = 5)

  expect_equal(chart2$ncomp, 2L)
  expect_equal(chart5$ncomp, 5L)
  expect_equal(length(chart2$eigenvalues), 2)
  expect_equal(length(chart5$eigenvalues), 5)
})

test_that("spm.phase1 stores .rust internals", {
  fd <- .test_spm_data()
  chart <- spm.phase1(fd, ncomp = 3)

  expect_true(!is.null(chart$.rust))
  expect_true(!is.null(chart$.rust$rotation))
  expect_true(!is.null(chart$.rust$scores))
  expect_true(!is.null(chart$.rust$mean))
})

test_that("spm.phase1 validates input", {
  expect_error(spm.phase1("not_fdata"), "fdata")
})

# =============================================================================
# spm.monitor
# =============================================================================

test_that("spm.monitor returns correct structure", {
  fd <- .test_spm_data(n = 60, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3)

  new_fd <- .test_spm_data(n = 20, seed = 99)
  result <- spm.monitor(chart, new_fd)

  expect_s3_class(result, "spm.monitor")
  expect_equal(length(result$t2), 20)
  expect_equal(length(result$spe), 20)
  expect_equal(length(result$t2.alarm), 20)
  expect_equal(length(result$spe.alarm), 20)
  expect_true(is.logical(result$t2.alarm))
  expect_true(is.logical(result$spe.alarm))
  expect_true(result$t2.ucl > 0)
  expect_true(result$spe.ucl > 0)
})

test_that("spm.monitor detects shifts", {
  fd <- .test_spm_data(n = 60, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3)

  # Create shifted data — large shift should trigger alarms
  set.seed(42)
  m <- 30
  argvals <- seq(0, 1, length.out = m)
  X_shift <- matrix(rnorm(20 * m, mean = 5), 20, m)
  fd_shift <- fdata(X_shift, argvals = argvals)
  result <- spm.monitor(chart, fd_shift)

  expect_true(any(result$t2.alarm) || any(result$spe.alarm))
})

test_that("spm.monitor validates input", {
  fd <- .test_spm_data()
  chart <- spm.phase1(fd, ncomp = 3)
  expect_error(spm.monitor("not_chart", fd), "spm.chart")
})

# =============================================================================
# spm.ewma
# =============================================================================

test_that("spm.ewma returns correct structure", {
  fd <- .test_spm_data(n = 60, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3)
  new_fd <- .test_spm_data(n = 20, seed = 99)

  result <- spm.ewma(chart, new_fd, lambda = 0.2)

  expect_true(!is.null(result$smoothed.scores))
  expect_true(!is.null(result$t2))
  expect_true(!is.null(result$spe))
  expect_equal(length(result$t2), 20)
  expect_equal(length(result$spe), 20)
  expect_true(is.logical(result$t2.alarm))
  expect_true(is.logical(result$spe.alarm))
  expect_true(result$t2.limit > 0)
  expect_true(result$spe.limit > 0)
})

test_that("spm.ewma lambda=1 gives raw scores", {
  fd <- .test_spm_data(n = 60, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3)
  new_fd <- .test_spm_data(n = 20, seed = 99)

  result_raw <- spm.monitor(chart, new_fd)
  result_ewma1 <- spm.ewma(chart, new_fd, lambda = 1.0)

  # SPE should be identical with lambda=1 (no smoothing on scores)
  expect_equal(result_ewma1$spe, result_raw$spe, tolerance = 1e-6)
})

test_that("spm.ewma validates input", {
  fd <- .test_spm_data()
  chart <- spm.phase1(fd, ncomp = 3)
  expect_error(spm.ewma("not_chart", fd), "spm.chart")
  expect_error(spm.ewma(chart, "not_fdata"), "fdata")
})

# =============================================================================
# spm.contributions
# =============================================================================

test_that("spm.contributions t2 returns matrix", {
  fd <- .test_spm_data(n = 60, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3)
  new_fd <- .test_spm_data(n = 20, seed = 99)
  result <- spm.monitor(chart, new_fd)

  m <- ncol(fd$data)
  contrib <- spm.contributions(result$scores, chart$eigenvalues,
                                grid.sizes = c(m), type = "t2")
  expect_true(is.matrix(contrib) || is.numeric(contrib))
  expect_true(nrow(contrib) == 20 || length(contrib) == 20)
})

test_that("spm.contributions validates eigenvalues", {
  expect_error(spm.contributions(matrix(1:6, 2, 3), type = "t2"),
               "eigenvalues")
})

# =============================================================================
# mfpca
# =============================================================================

test_that("mfpca returns correct structure", {
  set.seed(42)
  n <- 40
  fd1 <- fdata(matrix(rnorm(n * 20), n, 20), argvals = seq(0, 1, length.out = 20))
  fd2 <- fdata(matrix(rnorm(n * 15), n, 15), argvals = seq(0, 1, length.out = 15))

  result <- mfpca(list(fd1, fd2), ncomp = 3)

  expect_s3_class(result, "mfpca")
  expect_equal(nrow(result$scores), n)
  expect_equal(ncol(result$scores), 3)
  expect_equal(length(result$eigenvalues), 3)
  expect_true(all(result$eigenvalues > 0))
  expect_equal(length(result$eigenfunctions), 2)
  expect_equal(length(result$means), 2)
  expect_equal(result$grid.sizes, c(20, 15))
})

test_that("mfpca works with single variable", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(40 * 20), 40, 20), argvals = seq(0, 1, length.out = 20))
  result <- mfpca(list(fd), ncomp = 2)

  expect_s3_class(result, "mfpca")
  expect_equal(ncol(result$scores), 2)
})

test_that("mfpca validates input", {
  expect_error(mfpca("not_a_list"), "list")
  expect_error(mfpca(list()), "at least one")
  expect_error(mfpca(list("not_fdata")), "fdata")
})

# =============================================================================
# frcc.phase1 and frcc.monitor
# =============================================================================

test_that("frcc.phase1 returns correct structure", {
  set.seed(42)
  n <- 50
  m <- 20
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fd <- fdata(X, argvals = argvals)
  preds <- matrix(rnorm(n * 2), n, 2)

  chart <- frcc.phase1(fd, preds, ncomp = 3)

  expect_s3_class(chart, "frcc.chart")
  expect_true(chart$t2.ucl > 0)
  expect_true(chart$spe.ucl > 0)
  expect_equal(chart$ncomp, 3L)
})

test_that("frcc.monitor returns correct structure", {
  set.seed(42)
  n <- 50
  m <- 20
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fd <- fdata(X, argvals = argvals)
  preds <- matrix(rnorm(n * 2), n, 2)

  chart <- frcc.phase1(fd, preds, ncomp = 3)

  # New data
  X_new <- matrix(rnorm(10 * m), 10, m)
  fd_new <- fdata(X_new, argvals = argvals)
  preds_new <- matrix(rnorm(10 * 2), 10, 2)

  result <- frcc.monitor(chart, fd_new, preds_new)

  expect_s3_class(result, "spm.monitor")
  expect_equal(length(result$t2), 10)
  expect_equal(length(result$spe), 10)
})

# =============================================================================
# Print and plot methods
# =============================================================================

test_that("print.spm.chart works", {
  fd <- .test_spm_data()
  chart <- spm.phase1(fd, ncomp = 3)
  expect_output(print(chart), "SPM Control Chart")
})

test_that("print.spm.monitor works", {
  fd <- .test_spm_data(n = 60, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3)
  new_fd <- .test_spm_data(n = 20, seed = 99)
  result <- spm.monitor(chart, new_fd)
  expect_output(print(result), "SPM Monitoring Result")
})

test_that("print.mfpca works", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(40 * 20), 40, 20), argvals = seq(0, 1, length.out = 20))
  result <- mfpca(list(fd), ncomp = 2)
  expect_output(print(result), "Multivariate Functional PCA")
})

test_that("print.frcc.chart works", {
  set.seed(42)
  n <- 50; m <- 20
  fd <- fdata(matrix(rnorm(n * m), n, m), argvals = seq(0, 1, length.out = m))
  preds <- matrix(rnorm(n * 2), n, 2)
  chart <- frcc.phase1(fd, preds, ncomp = 3)
  expect_output(print(chart), "Functional Regression Control Chart")
})

test_that("plot.spm.chart runs without error", {
  fd <- .test_spm_data()
  chart <- spm.phase1(fd, ncomp = 3)
  pdf(NULL)  # suppress plot device
  expect_no_error(plot(chart))
  dev.off()
})

test_that("plot.spm.monitor runs without error", {
  fd <- .test_spm_data(n = 60, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3)
  new_fd <- .test_spm_data(n = 20, seed = 99)
  result <- spm.monitor(chart, new_fd)
  pdf(NULL)
  expect_no_error(plot(result))
  dev.off()
})
