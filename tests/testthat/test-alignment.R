test_that("srsf.transform returns correct dimensions", {
  fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
  q <- srsf.transform(fd)
  expect_s3_class(q, "fdata")
  expect_equal(nrow(q$data), 20)
  expect_equal(ncol(q$data), 10)
})

test_that("srsf round-trip approximately recovers original", {
  set.seed(42)
  argvals <- seq(0, 1, length.out = 50)
  f <- sin(2 * pi * argvals)
  fd <- fdata(matrix(f, nrow = 1), argvals = argvals)
  q <- srsf.transform(fd)
  f_hat <- srsf.inverse(q$data[1, ], argvals, f[1])
  # Approximate recovery (integration introduces small error)
  expect_true(cor(f, f_hat) > 0.95)
})

test_that("elastic.align returns correct structure", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(100), 10, 10), argvals = seq(0, 1, length.out = 10))
  res <- elastic.align(fd)
  expect_s3_class(res, "elastic.align")
  expect_s3_class(res$aligned, "fdata")
  expect_s3_class(res$gammas, "fdata")
  expect_equal(nrow(res$aligned$data), 10)
  expect_equal(length(res$distances), 10)
})

test_that("elastic.distance returns symmetric matrix with zero diagonal", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(60), 6, 10), argvals = seq(0, 1, length.out = 10))
  D <- elastic.distance(fd)
  expect_equal(nrow(D), 6)
  expect_equal(ncol(D), 6)
  # Symmetric
  expect_equal(D, t(D), tolerance = 1e-10)
  # Zero diagonal
  expect_equal(diag(D), rep(0, 6), tolerance = 1e-10)
})

test_that("elastic.distance cross-distance has correct dimensions", {
  set.seed(42)
  fd1 <- fdata(matrix(rnorm(60), 6, 10), argvals = seq(0, 1, length.out = 10))
  fd2 <- fdata(matrix(rnorm(40), 4, 10), argvals = seq(0, 1, length.out = 10))
  D <- elastic.distance(fd1, fd2)
  expect_equal(nrow(D), 6)
  expect_equal(ncol(D), 4)
})

test_that("karcher.mean returns correct structure", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(100), 10, 10), argvals = seq(0, 1, length.out = 10))
  km <- karcher.mean(fd, max.iter = 5)
  expect_s3_class(km, "karcher.mean")
  expect_s3_class(km$mean, "fdata")
  expect_equal(nrow(km$mean$data), 1)
  expect_equal(ncol(km$mean$data), 10)
  expect_true(km$n.iter >= 1)
  expect_type(km$converged, "logical")
})

test_that("metric dispatcher works with elastic", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(60), 6, 10), argvals = seq(0, 1, length.out = 10))
  D <- metric(fd, method = "elastic")
  expect_equal(nrow(D), 6)
  expect_equal(ncol(D), 6)
})

test_that("alignment input validation works", {
  expect_error(srsf.transform("not_fdata"))
  expect_error(elastic.align("not_fdata"))
  expect_error(elastic.distance("not_fdata"))
  expect_error(karcher.mean("not_fdata"))
})

test_that("print methods work without error", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(100), 10, 10), argvals = seq(0, 1, length.out = 10))
  res <- elastic.align(fd)
  expect_output(print(res), "Elastic Alignment")

  km <- karcher.mean(fd, max.iter = 3)
  expect_output(print(km), "Karcher Mean")
})

# =============================================================================
# Periodic alignment tests
# =============================================================================

test_that("periodic alignment improves VR on shifted periodic curves", {
  set.seed(42)
  m <- 100
  argvals <- seq(0, 2 * pi, length.out = m)
  n <- 15
  data <- matrix(0, n, m)
  for (i in 1:n) {
    shift <- sample(0:(m - 1), 1)
    idx <- ((seq_len(m) - 1L + shift) %% m) + 1L
    data[i, ] <- sin(argvals[idx])
  }
  fd <- fdata(data, argvals = argvals)

  # Without periodic: boundary constraints limit alignment
  res_std <- elastic.align(fd)
  var_orig <- mean(apply(fd$data, 2, var))
  vr_std <- 1 - mean(apply(res_std$aligned$data, 2, var)) / var_orig

  # With periodic: rotation + alignment should do much better
  res_per <- elastic.align(fd, periodic = TRUE)
  vr_per <- 1 - mean(apply(res_per$aligned$data, 2, var)) / var_orig

  expect_gt(vr_per, 0.90)
  expect_gt(vr_per, vr_std)
})

test_that("periodic.rotate returns correct structure and shifts", {
  set.seed(42)
  m <- 50
  argvals <- seq(0, 1, length.out = m)
  data <- matrix(0, 5, m)
  for (i in 1:5) {
    shift <- sample(0:(m - 1), 1)
    idx <- ((seq_len(m) - 1L + shift) %% m) + 1L
    data[i, ] <- sin(2 * pi * argvals[idx])
  }
  fd <- fdata(data, argvals = argvals)

  rot <- periodic.rotate(fd)
  expect_type(rot, "list")
  expect_s3_class(rot$fdataobj, "fdata")
  expect_type(rot$shifts, "integer")
  expect_equal(length(rot$shifts), 5)
  expect_equal(nrow(rot$fdataobj$data), 5)
  expect_equal(ncol(rot$fdataobj$data), m)
})

test_that("karcher.mean with periodic = TRUE converges", {
  set.seed(42)
  m <- 50
  argvals <- seq(0, 2 * pi, length.out = m)
  n <- 10
  data <- matrix(0, n, m)
  for (i in 1:n) {
    shift <- sample(0:(m - 1), 1)
    idx <- ((seq_len(m) - 1L + shift) %% m) + 1L
    data[i, ] <- sin(argvals[idx])
  }
  fd <- fdata(data, argvals = argvals)

  km <- karcher.mean(fd, max.iter = 10, periodic = TRUE)
  expect_s3_class(km, "karcher.mean")
  expect_false(is.null(km$rotations))
  expect_equal(length(km$rotations), n)
  expect_true(km$n.iter >= 1)
})

test_that("periodic = FALSE default leaves behavior unchanged", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(100), 10, 10), argvals = seq(0, 1, length.out = 10))

  res1 <- elastic.align(fd)
  res2 <- elastic.align(fd, periodic = FALSE)

  expect_equal(res1$aligned$data, res2$aligned$data)
  expect_null(res1$rotations)
  expect_null(res2$rotations)
})

test_that("rotations are integer grid shifts of correct length", {
  set.seed(42)
  m <- 80
  argvals <- seq(0, 2 * pi, length.out = m)
  n <- 8
  data <- matrix(0, n, m)
  for (i in 1:n) {
    shift <- sample(0:(m - 1), 1)
    idx <- ((seq_len(m) - 1L + shift) %% m) + 1L
    data[i, ] <- cos(argvals[idx])
  }
  fd <- fdata(data, argvals = argvals)

  res <- elastic.align(fd, periodic = TRUE)
  expect_type(res$rotations, "integer")
  expect_equal(length(res$rotations), n)

  # Print method should show Periodic: TRUE
  expect_output(print(res), "Periodic: TRUE")
})
