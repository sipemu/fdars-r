# Tests for functional mixed models

# =============================================================================
# fmm Tests
# =============================================================================

test_that("fmm returns correct structure", {
  set.seed(42)
  n_subjects <- 10
  n_obs_per <- 3
  n_total <- n_subjects * n_obs_per
  m <- 30
  t_grid <- seq(0, 1, length.out = m)

  X <- matrix(0, n_total, m)
  subject_ids <- rep(1:n_subjects, each = n_obs_per)
  for (i in 1:n_total) {
    subj <- subject_ids[i]
    X[i, ] <- sin(2 * pi * t_grid) + 0.3 * subj/n_subjects *
      cos(2 * pi * t_grid) + rnorm(m, sd = 0.1)
  }

  fd <- fdars::fdata(X, argvals = t_grid)
  result <- fdars::fmm(fd, subject_ids, ncomp = 2)

  expect_s3_class(result, "fmm")
  expect_true("mean.function" %in% names(result))
  expect_true("random.effects" %in% names(result))
  expect_true("fitted" %in% names(result))
  expect_true("residuals" %in% names(result))
  expect_true("sigma2.eps" %in% names(result))
  expect_true("n.subjects" %in% names(result))

  expect_s3_class(result$mean.function, "fdata")
  expect_s3_class(result$fitted, "fdata")
  expect_s3_class(result$residuals, "fdata")
  expect_equal(result$n.subjects, n_subjects)
})

test_that("fmm fitted values have correct dimensions", {
  set.seed(42)
  n_subjects <- 8
  n_obs_per <- 4
  n_total <- n_subjects * n_obs_per
  m <- 20
  t_grid <- seq(0, 1, length.out = m)

  X <- matrix(rnorm(n_total * m), n_total, m)
  subject_ids <- rep(1:n_subjects, each = n_obs_per)

  fd <- fdars::fdata(X, argvals = t_grid)
  result <- fdars::fmm(fd, subject_ids, ncomp = 2)

  expect_equal(nrow(result$fitted$data), n_total)
  expect_equal(ncol(result$fitted$data), m)
  expect_equal(nrow(result$residuals$data), n_total)
})

test_that("fmm with covariates works", {
  set.seed(42)
  n_subjects <- 10
  n_obs_per <- 3
  n_total <- n_subjects * n_obs_per
  m <- 20
  t_grid <- seq(0, 1, length.out = m)

  X <- matrix(rnorm(n_total * m), n_total, m)
  subject_ids <- rep(1:n_subjects, each = n_obs_per)
  covariates <- matrix(rnorm(n_total * 2), n_total, 2)

  fd <- fdars::fdata(X, argvals = t_grid)
  result <- fdars::fmm(fd, subject_ids, covariates = covariates, ncomp = 2)

  expect_s3_class(result, "fmm")
  expect_true("beta.functions" %in% names(result))
})

# =============================================================================
# fmm.predict Tests
# =============================================================================

test_that("fmm.predict returns fdata", {
  set.seed(42)
  n_subjects <- 8
  n_obs_per <- 3
  n_total <- n_subjects * n_obs_per
  m <- 20
  t_grid <- seq(0, 1, length.out = m)

  X <- matrix(rnorm(n_total * m), n_total, m)
  subject_ids <- rep(1:n_subjects, each = n_obs_per)
  covariates <- matrix(rnorm(n_total * 2), n_total, 2)

  fd <- fdars::fdata(X, argvals = t_grid)
  fit <- fdars::fmm(fd, subject_ids, covariates = covariates, ncomp = 2)

  new_cov <- matrix(rnorm(3 * 2), 3, 2)
  pred <- fdars::fmm.predict(fit, new.covariates = new_cov)

  expect_s3_class(pred, "fdata")
  expect_equal(nrow(pred$data), 3)
  expect_equal(ncol(pred$data), m)
})

test_that("fmm.predict without covariates returns mean", {
  set.seed(42)
  n_subjects <- 8
  n_obs_per <- 3
  n_total <- n_subjects * n_obs_per
  m <- 20
  t_grid <- seq(0, 1, length.out = m)

  X <- matrix(rnorm(n_total * m), n_total, m)
  subject_ids <- rep(1:n_subjects, each = n_obs_per)

  fd <- fdars::fdata(X, argvals = t_grid)
  fit <- fdars::fmm(fd, subject_ids, ncomp = 2)

  pred <- fdars::fmm.predict(fit)
  expect_s3_class(pred, "fdata")
})

# =============================================================================
# fmm.test.fixed Tests
# =============================================================================

test_that("fmm.test.fixed returns valid structure", {
  set.seed(42)
  n_subjects <- 10
  n_obs_per <- 3
  n_total <- n_subjects * n_obs_per
  m <- 20
  t_grid <- seq(0, 1, length.out = m)

  X <- matrix(rnorm(n_total * m), n_total, m)
  subject_ids <- rep(1:n_subjects, each = n_obs_per)
  covariates <- matrix(rnorm(n_total * 2), n_total, 2)

  fd <- fdars::fdata(X, argvals = t_grid)
  result <- fdars::fmm.test.fixed(fd, subject_ids, covariates,
                                    ncomp = 2, n.perm = 50, seed = 123)

  expect_s3_class(result, "fmm.test")
  expect_true("f.statistics" %in% names(result))
  expect_true("p.values" %in% names(result))
})

test_that("fmm.test.fixed p-values are valid", {
  set.seed(42)
  n_subjects <- 10
  n_obs_per <- 3
  n_total <- n_subjects * n_obs_per
  m <- 20
  t_grid <- seq(0, 1, length.out = m)

  X <- matrix(rnorm(n_total * m), n_total, m)
  subject_ids <- rep(1:n_subjects, each = n_obs_per)
  covariates <- matrix(rnorm(n_total), n_total, 1)

  fd <- fdars::fdata(X, argvals = t_grid)
  result <- fdars::fmm.test.fixed(fd, subject_ids, covariates,
                                    ncomp = 2, n.perm = 50, seed = 123)

  expect_true(all(result$p.values >= 0))
  expect_true(all(result$p.values <= 1))
})

# =============================================================================
# Input Validation Tests
# =============================================================================

test_that("fmm rejects non-fdata input", {
  X <- matrix(rnorm(100), 10, 10)
  expect_error(fdars::fmm(X, rep(1:5, each = 2)))
})

test_that("fmm rejects mismatched subject.ids length", {
  fd <- fdars::fdata(matrix(rnorm(100), 10, 10))
  expect_error(fdars::fmm(fd, 1:5), "Length of subject.ids")
})

# =============================================================================
# Print and Plot Tests
# =============================================================================

test_that("print.fmm produces output", {
  set.seed(42)
  n_subjects <- 6
  n_obs_per <- 3
  n_total <- n_subjects * n_obs_per
  m <- 15

  X <- matrix(rnorm(n_total * m), n_total, m)
  subject_ids <- rep(1:n_subjects, each = n_obs_per)

  fd <- fdars::fdata(X)
  fit <- fdars::fmm(fd, subject_ids, ncomp = 2)

  expect_output(print(fit), "Functional Mixed Model")
  expect_output(print(fit), "subjects")
})

test_that("print.fmm.test produces output", {
  set.seed(42)
  n_subjects <- 6
  n_obs_per <- 3
  n_total <- n_subjects * n_obs_per
  m <- 15

  X <- matrix(rnorm(n_total * m), n_total, m)
  subject_ids <- rep(1:n_subjects, each = n_obs_per)
  covariates <- matrix(rnorm(n_total), n_total, 1)

  fd <- fdars::fdata(X)
  result <- fdars::fmm.test.fixed(fd, subject_ids, covariates,
                                    ncomp = 2, n.perm = 20, seed = 42)

  expect_output(print(result), "Fixed Effects")
  expect_output(print(result), "Permutations")
})

test_that("plot.fmm works", {
  set.seed(42)
  n_subjects <- 6
  n_obs_per <- 3
  n_total <- n_subjects * n_obs_per
  m <- 15

  X <- matrix(rnorm(n_total * m), n_total, m)
  subject_ids <- rep(1:n_subjects, each = n_obs_per)

  fd <- fdars::fdata(X)
  fit <- fdars::fmm(fd, subject_ids, ncomp = 2)

  expect_no_error(plot(fit))
})
