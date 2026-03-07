# Tests for function-on-scalar regression

# =============================================================================
# fosr Tests
# =============================================================================

test_that("fosr returns correct structure", {
  set.seed(42)
  n <- 40
  m <- 30
  t_grid <- seq(0, 1, length.out = m)
  p <- 2

  X <- matrix(0, n, m)
  predictors <- matrix(rnorm(n * p), n, p)
  for (i in 1:n) {
    X[i, ] <- sin(2 * pi * t_grid) + predictors[i, 1] * t_grid +
      predictors[i, 2] * t_grid^2 + rnorm(m, sd = 0.1)
  }

  fd <- fdars::fdata(X, argvals = t_grid)
  result <- fdars::fosr(fd, predictors)

  expect_s3_class(result, "fosr")
  expect_true("intercept" %in% names(result))
  expect_true("beta" %in% names(result))
  expect_true("fitted" %in% names(result))
  expect_true("residuals" %in% names(result))
  expect_true("r.squared" %in% names(result))

  expect_s3_class(result$intercept, "fdata")
  expect_s3_class(result$beta, "fdata")
  expect_s3_class(result$fitted, "fdata")
  expect_s3_class(result$residuals, "fdata")
})

test_that("fosr R-squared is valid", {
  set.seed(42)
  n <- 40
  m <- 30
  t_grid <- seq(0, 1, length.out = m)

  X <- matrix(0, n, m)
  predictors <- matrix(rnorm(n), n, 1)
  for (i in 1:n) {
    X[i, ] <- predictors[i, 1] * sin(2 * pi * t_grid) + rnorm(m, sd = 0.1)
  }

  fd <- fdars::fdata(X, argvals = t_grid)
  result <- fdars::fosr(fd, predictors)

  expect_gte(result$r.squared, 0)
  expect_lte(result$r.squared, 1)
})

test_that("fosr with lambda regularization works", {
  set.seed(42)
  n <- 40
  m <- 30
  t_grid <- seq(0, 1, length.out = m)

  X <- matrix(rnorm(n * m), n, m)
  predictors <- matrix(rnorm(n * 2), n, 2)

  fd <- fdars::fdata(X, argvals = t_grid)

  result_noreg <- fdars::fosr(fd, predictors, lambda = 0)
  result_reg <- fdars::fosr(fd, predictors, lambda = 1)

  expect_s3_class(result_noreg, "fosr")
  expect_s3_class(result_reg, "fosr")
})

# =============================================================================
# fosr.fpc Tests
# =============================================================================

test_that("fosr.fpc returns correct structure", {
  set.seed(42)
  n <- 40
  m <- 30
  t_grid <- seq(0, 1, length.out = m)

  X <- matrix(rnorm(n * m), n, m)
  predictors <- matrix(rnorm(n * 2), n, 2)

  fd <- fdars::fdata(X, argvals = t_grid)
  result <- fdars::fosr.fpc(fd, predictors, ncomp = 3)

  expect_s3_class(result, "fosr")
  expect_true("ncomp" %in% names(result))
})

# =============================================================================
# fanova Tests
# =============================================================================

test_that("fanova returns correct structure", {
  set.seed(42)
  n_per_group <- 15
  m <- 30
  t_grid <- seq(0, 1, length.out = m)

  X <- matrix(0, 3 * n_per_group, m)
  for (i in 1:n_per_group) {
    X[i, ] <- sin(2 * pi * t_grid) + rnorm(m, sd = 0.1)
    X[n_per_group + i, ] <- cos(2 * pi * t_grid) + rnorm(m, sd = 0.1)
    X[2 * n_per_group + i, ] <- t_grid + rnorm(m, sd = 0.1)
  }
  groups <- rep(1:3, each = n_per_group)

  fd <- fdars::fdata(X, argvals = t_grid)
  result <- fdars::fanova(fd, groups, n.perm = 100)

  expect_s3_class(result, "fanova")
  expect_true("group.means" %in% names(result))
  expect_true("p.value" %in% names(result))
  expect_true("global.statistic" %in% names(result))
  expect_true("n.groups" %in% names(result))
})

test_that("fanova p-value is valid", {
  set.seed(42)
  n_per_group <- 15
  m <- 30
  t_grid <- seq(0, 1, length.out = m)

  X <- matrix(0, 3 * n_per_group, m)
  for (i in 1:n_per_group) {
    X[i, ] <- sin(2 * pi * t_grid) + rnorm(m, sd = 0.1)
    X[n_per_group + i, ] <- cos(2 * pi * t_grid) + rnorm(m, sd = 0.1)
    X[2 * n_per_group + i, ] <- t_grid + rnorm(m, sd = 0.1)
  }
  groups <- rep(1:3, each = n_per_group)

  fd <- fdars::fdata(X, argvals = t_grid)
  result <- fdars::fanova(fd, groups, n.perm = 100)

  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
})

test_that("fanova detects significant group differences", {
  set.seed(42)
  n_per_group <- 20
  m <- 30
  t_grid <- seq(0, 1, length.out = m)

  # Well-separated groups
  X <- matrix(0, 2 * n_per_group, m)
  for (i in 1:n_per_group) {
    X[i, ] <- sin(2 * pi * t_grid) + rnorm(m, sd = 0.1)
    X[n_per_group + i, ] <- 3 + cos(2 * pi * t_grid) + rnorm(m, sd = 0.1)
  }
  groups <- rep(1:2, each = n_per_group)

  fd <- fdars::fdata(X, argvals = t_grid)
  result <- fdars::fanova(fd, groups, n.perm = 200)

  expect_lt(result$p.value, 0.05)
})

# =============================================================================
# predict.fosr Tests
# =============================================================================

test_that("predict.fosr returns correct dimensions", {
  set.seed(42)
  n <- 40
  m <- 30
  t_grid <- seq(0, 1, length.out = m)

  X <- matrix(rnorm(n * m), n, m)
  predictors <- matrix(rnorm(n * 2), n, 2)

  fd <- fdars::fdata(X, argvals = t_grid)
  fit <- fdars::fosr(fd, predictors)

  # Predict with new data
  new_pred <- matrix(rnorm(5 * 2), 5, 2)
  pred <- predict(fit, new_pred)

  expect_s3_class(pred, "fdata")
  expect_equal(nrow(pred$data), 5)
  expect_equal(ncol(pred$data), m)
})

test_that("predict.fosr without new data returns fitted", {
  set.seed(42)
  n <- 30
  m <- 20
  X <- matrix(rnorm(n * m), n, m)
  predictors <- matrix(rnorm(n), n, 1)

  fd <- fdars::fdata(X)
  fit <- fdars::fosr(fd, predictors)

  pred <- predict(fit)
  expect_equal(pred$data, fit$fitted$data)
})

# =============================================================================
# Input Validation Tests
# =============================================================================

test_that("fosr rejects non-fdata input", {
  X <- matrix(rnorm(100), 10, 10)
  expect_error(fdars::fosr(X, matrix(rnorm(10), 10, 1)))
})

test_that("fosr rejects mismatched rows", {
  fd <- fdars::fdata(matrix(rnorm(100), 10, 10))
  expect_error(fdars::fosr(fd, matrix(rnorm(20), 20, 1)), "Number of rows")
})

test_that("fanova rejects mismatched groups length", {
  fd <- fdars::fdata(matrix(rnorm(100), 10, 10))
  expect_error(fdars::fanova(fd, 1:5), "Length of groups")
})

# =============================================================================
# Print and Plot Tests
# =============================================================================

test_that("print.fosr produces output", {
  set.seed(42)
  n <- 30
  m <- 20
  X <- matrix(rnorm(n * m), n, m)
  predictors <- matrix(rnorm(n), n, 1)

  fd <- fdars::fdata(X)
  fit <- fdars::fosr(fd, predictors)

  expect_output(print(fit), "Function-on-Scalar")
  expect_output(print(fit), "R-squared")
})

test_that("print.fanova produces output", {
  set.seed(42)
  n <- 30
  m <- 20
  X <- matrix(rnorm(n * m), n, m)
  groups <- rep(1:3, each = 10)

  fd <- fdars::fdata(X)
  result <- fdars::fanova(fd, groups, n.perm = 50)

  expect_output(print(result), "Functional ANOVA")
  expect_output(print(result), "P-value")
})

test_that("plot.fosr works", {
  set.seed(42)
  n <- 30
  m <- 20
  X <- matrix(rnorm(n * m), n, m)
  predictors <- matrix(rnorm(n * 2), n, 2)

  fd <- fdars::fdata(X)
  fit <- fdars::fosr(fd, predictors)

  expect_no_error(plot(fit))
})

test_that("plot.fanova works", {
  set.seed(42)
  n <- 30
  m <- 20
  X <- matrix(rnorm(n * m), n, m)
  groups <- rep(1:3, each = 10)

  fd <- fdars::fdata(X)
  result <- fdars::fanova(fd, groups, n.perm = 50)

  expect_no_error(plot(result))
})
