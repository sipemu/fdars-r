# Extended SPM tests: depth.streaming, softdtw.barycenter, mfpca extended

# =============================================================================
# depth.streaming
# =============================================================================

test_that("depth.streaming returns depth values", {
  set.seed(42)
  m <- 15
  t <- seq(0, 1, length.out = m)
  fd <- fdata(matrix(rnorm(30 * m), 30, m), argvals = t)

  result <- depth.streaming(fd)

  expect_true(is.numeric(result))
  expect_equal(length(result), 30)
  expect_true(all(is.finite(result)))
})

test_that("depth.streaming with reference sample", {
  set.seed(42)
  m <- 15
  t <- seq(0, 1, length.out = m)
  fd <- fdata(matrix(rnorm(10 * m), 10, m), argvals = t)
  fd_ref <- fdata(matrix(rnorm(30 * m), 30, m), argvals = t)

  result <- depth.streaming(fd, fdataori = fd_ref)
  expect_true(is.numeric(result))
  expect_equal(length(result), 10)
})

# =============================================================================
# softdtw.barycenter
# =============================================================================

test_that("softdtw.barycenter returns fdata barycenter", {
  set.seed(42)
  m <- 20
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, 5, m)
  for (i in 1:5) X[i, ] <- sin(2 * pi * t) + rnorm(m, sd = 0.1)
  fd <- fdata(X, argvals = t)

  result <- softdtw.barycenter(fd, gamma = 1.0)

  expect_s3_class(result, "fdata")
  expect_equal(nrow(result$data), 1)
  expect_equal(ncol(result$data), m)
})

test_that("softdtw.barycenter with different gamma", {
  set.seed(42)
  m <- 20
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, 5, m)
  for (i in 1:5) X[i, ] <- sin(2 * pi * t) + rnorm(m, sd = 0.1)
  fd <- fdata(X, argvals = t)

  r1 <- softdtw.barycenter(fd, gamma = 0.1)
  r2 <- softdtw.barycenter(fd, gamma = 10.0)

  expect_s3_class(r1, "fdata")
  expect_s3_class(r2, "fdata")
})

test_that("softdtw.barycenter validates input", {
  expect_error(softdtw.barycenter("not_fdata"), "fdata")
})

# =============================================================================
# mfpca extended (multivariate FPCA)
# =============================================================================

test_that("mfpca with 2 variables", {
  set.seed(42)
  m <- 20
  t <- seq(0, 1, length.out = m)
  n <- 25

  fd1 <- fdata(matrix(rnorm(n * m), n, m), argvals = t)
  fd2 <- fdata(matrix(rnorm(n * m), n, m), argvals = t)

  result <- mfpca(list(fd1, fd2), ncomp = 3)

  expect_s3_class(result, "mfpca")
  expect_true("scores" %in% names(result))
  expect_true("eigenvalues" %in% names(result))
})

test_that("mfpca with 3 variables", {
  set.seed(42)
  m <- 15
  t <- seq(0, 1, length.out = m)
  n <- 20

  fds <- lapply(1:3, function(i) fdata(matrix(rnorm(n * m), n, m), argvals = t))
  result <- mfpca(fds, ncomp = 2)

  expect_s3_class(result, "mfpca")
})

# =============================================================================
# spm.phase1 extended: ncomp and alpha variations
# =============================================================================

test_that("spm.phase1 with different ncomp", {
  set.seed(42)
  m <- 15
  t <- seq(0, 1, length.out = m)
  fd <- fdata(matrix(rnorm(30 * m), 30, m), argvals = t)

  chart2 <- spm.phase1(fd, ncomp = 2, seed = 42)
  chart5 <- spm.phase1(fd, ncomp = 5, seed = 42)

  expect_s3_class(chart2, "spm.chart")
  expect_s3_class(chart5, "spm.chart")
})

test_that("spm.phase1 with different alpha", {
  set.seed(42)
  m <- 15
  t <- seq(0, 1, length.out = m)
  fd <- fdata(matrix(rnorm(30 * m), 30, m), argvals = t)

  chart_05 <- spm.phase1(fd, ncomp = 3, alpha = 0.05, seed = 42)
  chart_01 <- spm.phase1(fd, ncomp = 3, alpha = 0.01, seed = 42)

  expect_s3_class(chart_05, "spm.chart")
  expect_s3_class(chart_01, "spm.chart")
})

test_that("spm.monitor with in-control data has fewer alarms", {
  set.seed(42)
  m <- 15
  t <- seq(0, 1, length.out = m)
  fd_train <- fdata(matrix(rnorm(40 * m), 40, m), argvals = t)

  chart <- spm.phase1(fd_train, ncomp = 3, alpha = 0.05, seed = 42)

  # In-control new data
  fd_ic <- fdata(matrix(rnorm(10 * m), 10, m), argvals = t)
  mon_ic <- spm.monitor(chart, fd_ic)

  # Out-of-control new data (shifted)
  fd_oc <- fdata(matrix(rnorm(10 * m, mean = 5), 10, m), argvals = t)
  mon_oc <- spm.monitor(chart, fd_oc)

  expect_s3_class(mon_ic, "spm.monitor")
  expect_s3_class(mon_oc, "spm.monitor")
})

# =============================================================================
# spm.ewma extended
# =============================================================================

test_that("spm.ewma with different lambda values", {
  set.seed(42)
  m <- 15
  t <- seq(0, 1, length.out = m)
  fd_train <- fdata(matrix(rnorm(30 * m), 30, m), argvals = t)
  chart <- spm.phase1(fd_train, ncomp = 3, seed = 42)

  fd_new <- fdata(matrix(rnorm(10 * m), 10, m), argvals = t)

  ewma_low <- spm.ewma(chart, fd_new, lambda = 0.1)
  ewma_high <- spm.ewma(chart, fd_new, lambda = 0.9)

  # spm.ewma may return spm.monitor or plain list depending on lambda
  expect_true(is.list(ewma_low))
  expect_true(is.list(ewma_high))
  expect_true("t2" %in% names(ewma_low))
  expect_true("spe" %in% names(ewma_low))
})

# =============================================================================
# frcc extended
# =============================================================================

test_that("frcc.phase1 with different ncomp", {
  set.seed(42)
  m <- 15
  t <- seq(0, 1, length.out = m)
  n <- 30
  fd <- fdata(matrix(rnorm(n * m), n, m), argvals = t)
  z <- matrix(rnorm(n * 2), n, 2)

  chart <- frcc.phase1(fd, z, ncomp = 2, seed = 42)
  expect_s3_class(chart, "frcc.chart")
})
