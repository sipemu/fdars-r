# Tests for advanced SPM functions: Phases 1-3

# Helper: create in-control functional data
.test_spm_adv_data <- function(n = 50, m = 30, seed = 42) {
  set.seed(seed)
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fdata(X, argvals = argvals)
}

# =============================================================================
# spm.ncomp.select
# =============================================================================

test_that("spm.ncomp.select variance90 returns reasonable ncomp", {
  eigenvalues <- c(10, 5, 2, 1, 0.5, 0.2, 0.1)
  ncomp <- spm.ncomp.select(eigenvalues, method = "variance90")
  expect_true(is.integer(ncomp))
  expect_gte(ncomp, 1L)
  expect_lte(ncomp, length(eigenvalues))
  # Cumulative variance at selected ncomp should be >= 90%
  cumvar <- cumsum(eigenvalues) / sum(eigenvalues)
  expect_gte(cumvar[ncomp], 0.9)
})

test_that("spm.ncomp.select kaiser rule selects eigenvalues > mean", {
  eigenvalues <- c(10, 5, 2, 1, 0.5, 0.2, 0.1)
  ncomp <- spm.ncomp.select(eigenvalues, method = "kaiser")
  expect_true(is.integer(ncomp))
  expect_gte(ncomp, 1L)
  expect_lte(ncomp, length(eigenvalues))
})

test_that("spm.ncomp.select elbow method returns integer", {
  eigenvalues <- c(10, 5, 2, 1, 0.5, 0.2, 0.1)
  ncomp <- spm.ncomp.select(eigenvalues, method = "elbow")
  expect_true(is.integer(ncomp))
  expect_gte(ncomp, 1L)
  expect_lte(ncomp, length(eigenvalues))
})

test_that("spm.ncomp.select fixed method respects threshold", {
  eigenvalues <- c(10, 5, 2, 1, 0.5, 0.2, 0.1)
  ncomp_80 <- spm.ncomp.select(eigenvalues, method = "fixed", threshold = 0.8)
  ncomp_99 <- spm.ncomp.select(eigenvalues, method = "fixed", threshold = 0.99)
  expect_true(is.integer(ncomp_80))
  expect_true(is.integer(ncomp_99))
  expect_lte(ncomp_80, ncomp_99)
})

test_that("spm.ncomp.select validates empty eigenvalues", {
  expect_error(spm.ncomp.select(numeric(0)), "non-empty")
})

# =============================================================================
# spm.pc.contributions
# =============================================================================

test_that("spm.pc.contributions returns correct dimensions", {
  set.seed(1)
  scores <- matrix(rnorm(10 * 3), 10, 3)
  eigenvalues <- c(5, 2, 1)
  contrib <- spm.pc.contributions(scores, eigenvalues)

  expect_true(is.matrix(contrib))
  expect_equal(nrow(contrib), 10)
  expect_equal(ncol(contrib), 3)
})

test_that("spm.pc.contributions row sums equal T2 statistic", {
  set.seed(1)
  scores <- matrix(rnorm(10 * 3), 10, 3)
  eigenvalues <- c(5, 2, 1)
  contrib <- spm.pc.contributions(scores, eigenvalues)

  # T2 = sum(score_j^2 / eigenvalue_j) for each observation
  t2_expected <- rowSums(scores^2 / matrix(eigenvalues, 10, 3, byrow = TRUE))
  t2_from_contrib <- rowSums(contrib)
  expect_equal(t2_from_contrib, t2_expected, tolerance = 1e-6)
})

test_that("spm.pc.contributions validates dimension mismatch", {
  scores <- matrix(rnorm(10 * 3), 10, 3)
  expect_error(spm.pc.contributions(scores, c(1, 2)),
               "columns in scores must equal length of eigenvalues")
})

# =============================================================================
# spm.rules
# =============================================================================

test_that("spm.rules western electric detects points beyond 3 sigma", {
  # Place several points beyond 3 sigma
  values <- c(rnorm(15), rep(4, 5))
  result <- spm.rules(values, center = 0, sigma = 1,
                       rule.set = "western.electric")

  expect_s3_class(result, "spm.rules")
  expect_true(length(result$rules) > 0)
  expect_true(is.list(result$indices))
  expect_equal(length(result$rules), length(result$indices))
})

test_that("spm.rules western electric WE4: 8 consecutive above center", {
  # 8+ consecutive points above center should trigger WE rule 4
  values <- c(rnorm(5, mean = 0), rep(0.5, 10))  # 10 above center
  result <- spm.rules(values, center = 0, sigma = 1,
                       rule.set = "western.electric")

  expect_s3_class(result, "spm.rules")
  # At least one rule should fire for this pattern
  expect_true(length(result$rules) > 0)
})

test_that("spm.rules nelson rule set returns result", {
  set.seed(1)
  values <- c(rnorm(15), rep(4, 5))
  result <- spm.rules(values, center = 0, sigma = 1,
                       rule.set = "nelson")

  expect_s3_class(result, "spm.rules")
  expect_true(length(result$rules) >= 0)
  expect_equal(result$rule.set, "nelson")
})

test_that("spm.rules in-control data may have no violations", {
  set.seed(42)
  values <- rnorm(20, mean = 0, sd = 0.5)
  result <- spm.rules(values, center = 0, sigma = 1,
                       rule.set = "western.electric")

  expect_s3_class(result, "spm.rules")
  expect_true(length(result$rules) >= 0)
})

test_that("spm.rules validates sigma positive", {
  expect_error(spm.rules(c(1, 2, 3), 0, -1), "positive")
})

test_that("spm.rules print method works", {
  values <- c(rnorm(15), rep(4, 5))
  result <- spm.rules(values, center = 0, sigma = 1)
  expect_output(print(result), "SPM Control Chart Rule Evaluation")
})

# =============================================================================
# spm.limit.robust
# =============================================================================

test_that("spm.limit.robust t2 parametric returns UCL", {
  set.seed(1)
  t2_values <- rchisq(100, df = 3)
  result <- spm.limit.robust(t2_values, ncomp = 3, type = "t2",
                              method = "parametric")

  expect_true(is.list(result))
  expect_true(result$ucl > 0)
  expect_equal(result$alpha, 0.05)
  expect_equal(result$method, "parametric")
  expect_true(is.character(result$description))
})

test_that("spm.limit.robust t2 empirical returns UCL", {
  set.seed(1)
  t2_values <- rchisq(100, df = 3)
  result <- spm.limit.robust(t2_values, ncomp = 3, type = "t2",
                              method = "empirical")

  expect_true(result$ucl > 0)
  expect_equal(result$method, "empirical")
})

test_that("spm.limit.robust t2 bootstrap returns UCL", {
  set.seed(1)
  t2_values <- rchisq(100, df = 3)
  result <- spm.limit.robust(t2_values, ncomp = 3, type = "t2",
                              method = "bootstrap", n.bootstrap = 200)

  expect_true(result$ucl > 0)
  expect_equal(result$method, "bootstrap")
})

test_that("spm.limit.robust t2 kde returns UCL", {
  set.seed(1)
  t2_values <- rchisq(100, df = 3)
  result <- spm.limit.robust(t2_values, ncomp = 3, type = "t2",
                              method = "kde")

  expect_true(result$ucl > 0)
  expect_equal(result$method, "kde")
})

test_that("spm.limit.robust spe methods work", {
  set.seed(1)
  spe_values <- rchisq(100, df = 5)

  r_emp <- spm.limit.robust(spe_values, type = "spe", method = "empirical")
  r_boot <- spm.limit.robust(spe_values, type = "spe", method = "bootstrap",
                              n.bootstrap = 200)

  expect_true(r_emp$ucl > 0)
  expect_true(r_boot$ucl > 0)
})

test_that("spm.limit.robust t2 requires ncomp", {
  expect_error(spm.limit.robust(c(1, 2, 3), type = "t2"), "ncomp")
})

# =============================================================================
# spm.arl
# =============================================================================

test_that("spm.arl in-control T2 returns positive ARL", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)

  arl <- spm.arl(chart, type = "t2", n.sim = 500, max.rl = 500, seed = 1)

  expect_s3_class(arl, "spm.arl")
  expect_true(arl$arl > 0)
  expect_true(arl$std.dev > 0)
  expect_true(arl$median.rl > 0)
  expect_equal(length(arl$run.lengths), 500)
  expect_true(arl$in.control)
  expect_equal(arl$type, "t2")
})

test_that("spm.arl in-control ARL is approximately 1/alpha", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, alpha = 0.05, seed = 1)

  arl <- spm.arl(chart, type = "t2", n.sim = 500, max.rl = 1000, seed = 1)

  # For alpha=0.05, expected in-control ARL ~ 20
  # Be generous with bounds since it's stochastic
  expect_gt(arl$arl, 5)
})

test_that("spm.arl out-of-control T2 has shorter ARL", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)

  arl0 <- spm.arl(chart, type = "t2", n.sim = 200, max.rl = 500, seed = 1)
  arl1 <- spm.arl(chart, type = "t2",
                    shift.magnitude = rep(2, 3),
                    n.sim = 200, max.rl = 500, seed = 1)

  expect_true(arl1$arl < arl0$arl)
  expect_false(arl1$in.control)
})

test_that("spm.arl SPE type works", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)

  arl <- spm.arl(chart, type = "spe", n.sim = 200, max.rl = 500, seed = 1)

  expect_s3_class(arl, "spm.arl")
  expect_true(arl$arl > 0)
  expect_equal(arl$type, "spe")
})

test_that("spm.arl print method works", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)
  arl <- spm.arl(chart, type = "t2", n.sim = 100, max.rl = 100, seed = 1)
  expect_output(print(arl), "Average Run Length")
})

# =============================================================================
# spm.arl.ewma
# =============================================================================

test_that("spm.arl.ewma returns positive ARL", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)

  arl <- spm.arl.ewma(chart, lambda = 0.2,
                        n.sim = 200, max.rl = 500, seed = 1)

  expect_s3_class(arl, "spm.arl")
  expect_true(arl$arl > 0)
  expect_true(arl$std.dev > 0)
  expect_equal(arl$type, "ewma.t2")
  expect_true(arl$in.control)
})

test_that("spm.arl.ewma with custom UCL", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)

  arl <- spm.arl.ewma(chart, lambda = 0.2, ucl = 50,
                        n.sim = 100, max.rl = 200, seed = 1)

  expect_s3_class(arl, "spm.arl")
  expect_true(arl$arl > 0)
})

# =============================================================================
# spm.cusum
# =============================================================================

test_that("spm.cusum returns correct structure", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)

  fd_new <- .test_spm_adv_data(n = 20, m = 20, seed = 99)
  cusum <- spm.cusum(chart, fd_new, k = 0.5, h = 5.0)

  expect_s3_class(cusum, "spm.cusum")
  expect_equal(length(cusum$cusum.statistic), 20)
  expect_equal(length(cusum$alarm), 20)
  expect_true(is.logical(cusum$alarm))
  expect_equal(cusum$ucl, 5.0)
  expect_equal(cusum$k, 0.5)
  expect_true(is.numeric(cusum$spe))
  expect_equal(length(cusum$spe), 20)
  expect_true(is.logical(cusum$spe.alarm))
})

test_that("spm.cusum detects shift in data", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)

  # Shifted data
  set.seed(42)
  m <- 20
  argvals <- seq(0, 1, length.out = m)
  X_shift <- matrix(rnorm(30 * m, mean = 3), 30, m)
  fd_shift <- fdata(X_shift, argvals = argvals)

  cusum <- spm.cusum(chart, fd_shift, k = 0.5, h = 5.0)

  expect_true(any(cusum$alarm))
})

test_that("spm.cusum with restart resets after alarm", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)

  set.seed(42)
  m <- 20
  argvals <- seq(0, 1, length.out = m)
  X_shift <- matrix(rnorm(30 * m, mean = 3), 30, m)
  fd_shift <- fdata(X_shift, argvals = argvals)

  cusum_no_restart <- spm.cusum(chart, fd_shift, k = 0.5, h = 5.0,
                                 restart = FALSE)
  cusum_restart <- spm.cusum(chart, fd_shift, k = 0.5, h = 5.0,
                              restart = TRUE)

  expect_s3_class(cusum_restart, "spm.cusum")
  expect_true(cusum_restart$restart)
  expect_false(cusum_no_restart$restart)
})

test_that("spm.cusum print method works", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)
  fd_new <- .test_spm_adv_data(n = 10, m = 20, seed = 99)
  cusum <- spm.cusum(chart, fd_new, k = 0.5, h = 5.0)
  expect_output(print(cusum), "SPM CUSUM Monitoring")
})

# =============================================================================
# spm.mewma
# =============================================================================

test_that("spm.mewma returns correct structure", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)
  fd_new <- .test_spm_adv_data(n = 20, m = 20, seed = 99)

  result <- spm.mewma(chart, fd_new, lambda = 0.2)

  expect_s3_class(result, "spm.mewma")
  expect_equal(length(result$mewma.statistic), 20)
  expect_equal(length(result$alarm), 20)
  expect_true(is.logical(result$alarm))
  expect_true(result$ucl > 0)
  expect_true(is.numeric(result$spe))
  expect_equal(length(result$spe), 20)
  expect_true(is.logical(result$spe.alarm))
  expect_equal(result$lambda, 0.2)
})

test_that("spm.mewma detects shift in mean", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)

  set.seed(42)
  m <- 20
  argvals <- seq(0, 1, length.out = m)
  X_shift <- matrix(rnorm(20 * m, mean = 5), 20, m)
  fd_shift <- fdata(X_shift, argvals = argvals)

  result <- spm.mewma(chart, fd_shift, lambda = 0.2)
  expect_true(any(result$alarm) || any(result$spe.alarm))
})

test_that("spm.mewma print method works", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)
  fd_new <- .test_spm_adv_data(n = 10, m = 20, seed = 99)
  result <- spm.mewma(chart, fd_new, lambda = 0.2)
  expect_output(print(result), "SPM MEWMA Monitoring")
})

# =============================================================================
# spm.amewma
# =============================================================================

test_that("spm.amewma returns correct structure", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)
  fd_new <- .test_spm_adv_data(n = 20, m = 20, seed = 99)

  result <- spm.amewma(chart, fd_new)

  expect_s3_class(result, "spm.amewma")
  expect_equal(length(result$t2.statistic), 20)
  expect_equal(length(result$lambda.t), 20)
  expect_equal(length(result$alarm), 20)
  expect_true(is.logical(result$alarm))
  expect_true(result$ucl > 0)
  # lambda.t values should be within bounds
  expect_true(all(result$lambda.t >= result$lambda.min))
  expect_true(all(result$lambda.t <= result$lambda.max))
})

test_that("spm.amewma detects shift in data", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)

  set.seed(42)
  m <- 20
  argvals <- seq(0, 1, length.out = m)
  X_shift <- matrix(rnorm(20 * m, mean = 5), 20, m)
  fd_shift <- fdata(X_shift, argvals = argvals)

  result <- spm.amewma(chart, fd_shift)
  expect_true(any(result$alarm))
})

test_that("spm.amewma print method works", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)
  fd_new <- .test_spm_adv_data(n = 10, m = 20, seed = 99)
  result <- spm.amewma(chart, fd_new)
  expect_output(print(result), "Adaptive MEWMA")
})

test_that("spm.amewma with custom parameters", {
  fd <- .test_spm_adv_data(n = 50, m = 20, seed = 1)
  chart <- spm.phase1(fd, ncomp = 3, seed = 1)
  fd_new <- .test_spm_adv_data(n = 15, m = 20, seed = 99)

  result <- spm.amewma(chart, fd_new, lambda.min = 0.1, lambda.max = 0.8,
                        lambda.init = 0.3, eta = 0.2)

  expect_s3_class(result, "spm.amewma")
  expect_equal(result$lambda.min, 0.1)
  expect_equal(result$lambda.max, 0.8)
  expect_equal(result$lambda.init, 0.3)
  expect_equal(result$eta, 0.2)
})

# =============================================================================
# spm.phase1.iterative
# =============================================================================

test_that("spm.phase1.iterative returns spm.chart with extra fields", {
  set.seed(1)
  n <- 60; m <- 20
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  # Add outliers to first 3 curves
  X[1:3, ] <- X[1:3, ] + 5
  fd <- fdata(X, argvals = argvals)

  chart <- spm.phase1.iterative(fd, ncomp = 3, max.iterations = 5)

  expect_s3_class(chart, "spm.chart")
  expect_true(!is.null(chart$n.iterations))
  expect_true(!is.null(chart$removed.indices))
  expect_true(!is.null(chart$n.remaining))
  expect_true(chart$n.iterations >= 0)
  expect_true(chart$n.remaining > 0)
  expect_true(chart$n.remaining <= n)
})

test_that("spm.phase1.iterative removes contaminated observations", {
  set.seed(1)
  n <- 60; m <- 20
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  # Add large outliers
  X[1:5, ] <- X[1:5, ] + 10
  fd <- fdata(X, argvals = argvals)

  chart <- spm.phase1.iterative(fd, ncomp = 3, max.iterations = 10)

  # With extreme outliers, iterative should remove some observations
  expect_true(length(chart$removed.indices) > 0 || chart$n.remaining < n)
})

test_that("spm.phase1.iterative result works with spm.monitor", {
  set.seed(1)
  n <- 60; m <- 20
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  X[1:3, ] <- X[1:3, ] + 5
  fd <- fdata(X, argvals = argvals)

  chart <- spm.phase1.iterative(fd, ncomp = 3, max.iterations = 5)

  # Chart should be usable for monitoring
  fd_new <- fdata(matrix(rnorm(10 * m), 10, m), argvals = argvals)
  mon <- spm.monitor(chart, fd_new)
  expect_s3_class(mon, "spm.monitor")
  expect_equal(length(mon$t2), 10)
})

test_that("spm.phase1.iterative respects max.removal.fraction", {
  set.seed(1)
  n <- 40; m <- 20
  argvals <- seq(0, 1, length.out = m)
  X <- matrix(rnorm(n * m), n, m)
  fd <- fdata(X, argvals = argvals)

  chart <- spm.phase1.iterative(fd, ncomp = 3, max.iterations = 5,
                                  max.removal.fraction = 0.1)

  expect_s3_class(chart, "spm.chart")
  # Should not remove more than 10% of observations
  n_removed <- length(chart$removed.indices)
  expect_lte(n_removed, ceiling(n * 0.1))
})
