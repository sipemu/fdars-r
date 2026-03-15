# Smoke tests for plot methods that lack test coverage

# =============================================================================
# Elastic FPCA plot methods
# =============================================================================

test_that("plot.vert.fpca works for all types", {
  set.seed(42)
  t_grid <- seq(0, 1, length.out = 20)
  n <- 20
  X <- matrix(0, n, 20)
  for (i in 1:n) X[i, ] <- sin(2 * pi * t_grid) + rnorm(20, sd = 0.1) * (1 + 0.3 * i / n)
  fd <- fdata(X, argvals = t_grid)
  km <- karcher.mean(fd, max.iter = 3)
  vfpca <- vert.fpca(km, ncomp = 2)

  pdf(NULL)
  expect_no_error(plot(vfpca, type = "scores"))
  expect_no_error(plot(vfpca, type = "eigenfunctions"))
  expect_no_error(plot(vfpca, type = "variance"))
  dev.off()
})

test_that("plot.horiz.fpca works for all types", {
  set.seed(42)
  t_grid <- seq(0, 1, length.out = 20)
  n <- 20
  X <- matrix(0, n, 20)
  for (i in 1:n) {
    phase <- rnorm(1, 0, 0.05)
    X[i, ] <- sin(2 * pi * (t_grid + phase)) + rnorm(20, sd = 0.05)
  }
  fd <- fdata(X, argvals = t_grid)
  km <- karcher.mean(fd, max.iter = 3)
  hfpca <- horiz.fpca(km, ncomp = 2)

  pdf(NULL)
  expect_no_error(plot(hfpca, type = "scores"))
  expect_no_error(plot(hfpca, type = "eigenfunctions"))
  expect_no_error(plot(hfpca, type = "warps"))
  expect_no_error(plot(hfpca, type = "variance"))
  dev.off()
})

test_that("plot.joint.fpca works", {
  set.seed(42)
  t_grid <- seq(0, 1, length.out = 20)
  n <- 20
  X <- matrix(0, n, 20)
  for (i in 1:n) {
    amp <- rnorm(1, 1, 0.2)
    phase <- rnorm(1, 0, 0.05)
    X[i, ] <- amp * sin(2 * pi * (t_grid + phase)) + rnorm(20, sd = 0.05)
  }
  fd <- fdata(X, argvals = t_grid)
  km <- karcher.mean(fd, max.iter = 3)
  jfpca <- joint.fpca(km, ncomp = 2)

  pdf(NULL)
  expect_no_error(plot(jfpca, type = "scores"))
  expect_no_error(plot(jfpca, type = "variance"))
  dev.off()
})

# =============================================================================
# Elastic logistic plot
# =============================================================================

test_that("plot.elastic.logistic works", {
  set.seed(42)
  t_grid <- seq(0, 1, length.out = 20)
  n <- 30
  X <- matrix(0, n, 20)
  for (i in 1:n) {
    amp <- if (i <= 15) rnorm(1, 1, 0.2) else rnorm(1, -1, 0.2)
    X[i, ] <- amp * sin(2 * pi * t_grid) + rnorm(20, sd = 0.1)
  }
  fd <- fdata(X, argvals = t_grid)
  y <- c(rep(0, 15), rep(1, 15))

  fit <- elastic.logistic(fd, y, max.iter = 3)

  pdf(NULL)
  expect_no_error(plot(fit))
  dev.off()
})

# =============================================================================
# Karcher mean plot
# =============================================================================

test_that("plot.karcher.mean works", {
  set.seed(42)
  t_grid <- seq(0, 1, length.out = 20)
  n <- 15
  X <- matrix(0, n, 20)
  for (i in 1:n) X[i, ] <- sin(2 * pi * t_grid) + rnorm(20, sd = 0.1)
  fd <- fdata(X, argvals = t_grid)
  km <- karcher.mean(fd, max.iter = 3)

  pdf(NULL)
  expect_no_error(plot(km))
  dev.off()
})

# =============================================================================
# Seasonal analysis plot methods
# =============================================================================

test_that("plot.lomb_scargle_result works", {
  set.seed(42)
  t_grid <- sort(runif(80, 0, 10))
  y <- sin(2 * pi * t_grid / 2.5) + rnorm(80, sd = 0.3)
  fd <- fdata(matrix(y, nrow = 1), argvals = t_grid)
  result <- lomb.scargle(fd)

  pdf(NULL)
  expect_no_error(plot(result))
  dev.off()
})

test_that("plot.matrix_profile_result works", {
  set.seed(42)
  n <- 100
  t_grid <- seq(0, 1, length.out = n)
  x <- sin(seq(0, 8 * pi, length.out = n)) + rnorm(n, sd = 0.1)
  fd <- fdata(matrix(x, nrow = 1), argvals = t_grid)
  result <- matrix.profile(fd, subsequence_length = 10)

  pdf(NULL)
  expect_no_error(plot(result))
  dev.off()
})

test_that("plot.peak_detection works", {
  set.seed(42)
  t_grid <- seq(0, 1, length.out = 100)
  x <- sin(seq(0, 4 * pi, length.out = 100)) + rnorm(100, sd = 0.1)
  fd <- fdata(matrix(x, nrow = 1), argvals = t_grid)
  result <- detect.peaks(fd)

  pdf(NULL)
  expect_no_error(plot(result))
  dev.off()
})

# =============================================================================
# FOSR plot
# =============================================================================

test_that("plot.fosr works", {
  set.seed(42)
  n <- 30; m <- 15
  t_grid <- seq(0, 1, length.out = m)
  Y <- matrix(rnorm(n * m), n, m)
  fd <- fdata(Y, argvals = t_grid)
  X <- cbind(rnorm(n), rnorm(n))

  fit <- fosr(fd, X)
  pdf(NULL)
  expect_no_error(plot(fit))
  dev.off()
})

# =============================================================================
# Alignment quality plot
# =============================================================================

test_that("plot.alignment.quality works", {
  set.seed(42)
  t_grid <- seq(0, 1, length.out = 20)
  n <- 15
  X <- matrix(0, n, 20)
  for (i in 1:n) {
    phase <- rnorm(1, 0, 0.05)
    X[i, ] <- sin(2 * pi * (t_grid + phase)) + rnorm(20, sd = 0.05)
  }
  fd <- fdata(X, argvals = t_grid)
  km <- karcher.mean(fd, max.iter = 3)
  aq <- alignment.quality(fd, km)

  pdf(NULL)
  expect_no_error(plot(aq))
  dev.off()
})

# =============================================================================
# fregre.bootstrap.ci plot
# =============================================================================

test_that("plot.fregre.bootstrap.ci works", {
  set.seed(42)
  t_grid <- seq(0, 1, length.out = 15)
  n <- 25
  X <- matrix(0, n, 15)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1, 1, 0.3)
    X[i, ] <- amp * sin(2 * pi * t_grid) + rnorm(15, sd = 0.1)
    y[i] <- 2 * amp + rnorm(1, sd = 0.2)
  }
  fd <- fdata(X, argvals = t_grid)
  model <- fregre.lm(fd, y, ncomp = 3)
  ci <- fregre.bootstrap.ci(model, n.boot = 20, seed = 42)

  pdf(NULL)
  expect_no_error(plot(ci))
  dev.off()
})
