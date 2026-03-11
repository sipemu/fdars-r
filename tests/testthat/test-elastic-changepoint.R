# Tests for elastic changepoint detection

test_that("elastic.changepoint amplitude returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 20
  X <- matrix(0, n, 20)
  # First half: sin, second half: cos (amplitude change)
  for (i in 1:10) X[i, ] <- sin(2 * pi * t) + rnorm(20, sd = 0.1)
  for (i in 11:20) X[i, ] <- 2 * sin(2 * pi * t) + rnorm(20, sd = 0.1)
  fd <- fdata(X, argvals = t)

  result <- elastic.changepoint(fd, type = "amplitude", max.iter = 3,
                                 n.mc = 50, seed = 42)

  expect_s3_class(result, "elastic.changepoint")
  expect_true(is.finite(result$changepoint))
  expect_true(result$changepoint >= 1 && result$changepoint <= n)
  expect_true(is.finite(result$test.statistic))
  expect_true(result$test.statistic >= 0)
  expect_true(is.finite(result$p.value))
  expect_gte(result$p.value, 0)
  expect_lte(result$p.value, 1)
  expect_equal(result$type, "amplitude")
  expect_equal(result$n, n)
  expect_true(length(result$cusum.values) > 0)
})

test_that("elastic.changepoint phase returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 20
  X <- matrix(0, n, 20)
  for (i in 1:10) X[i, ] <- sin(2 * pi * t) + rnorm(20, sd = 0.1)
  for (i in 11:20) X[i, ] <- sin(2 * pi * (t + 0.2)) + rnorm(20, sd = 0.1)
  fd <- fdata(X, argvals = t)

  result <- elastic.changepoint(fd, type = "phase", max.iter = 3,
                                 n.mc = 50, seed = 42)

  expect_s3_class(result, "elastic.changepoint")
  expect_equal(result$type, "phase")
  expect_true(is.finite(result$changepoint))
  expect_true(is.finite(result$p.value))
})

test_that("elastic.changepoint fpca returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 20
  X <- matrix(0, n, 20)
  for (i in 1:10) X[i, ] <- sin(2 * pi * t) + rnorm(20, sd = 0.1)
  for (i in 11:20) X[i, ] <- cos(2 * pi * t) + rnorm(20, sd = 0.1)
  fd <- fdata(X, argvals = t)

  result <- elastic.changepoint(fd, type = "fpca", pca.method = "vertical",
                                 ncomp = 2, max.iter = 3, n.mc = 50, seed = 42)

  expect_s3_class(result, "elastic.changepoint")
  expect_equal(result$type, "fpca")
  expect_true(is.finite(result$changepoint))
  expect_true(is.finite(result$p.value))
})

test_that("elastic.changepoint input validation works", {
  expect_error(elastic.changepoint("not_fdata"))
})

test_that("print.elastic.changepoint works", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 20
  X <- matrix(0, n, 20)
  for (i in 1:10) X[i, ] <- sin(2 * pi * t) + rnorm(20, sd = 0.1)
  for (i in 11:20) X[i, ] <- 2 * sin(2 * pi * t) + rnorm(20, sd = 0.1)
  fd <- fdata(X, argvals = t)

  result <- elastic.changepoint(fd, type = "amplitude", max.iter = 3,
                                 n.mc = 50, seed = 42)

  expect_output(print(result), "Elastic Changepoint Detection")
  expect_output(print(result), "amplitude")
})

test_that("plot.elastic.changepoint works", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 20
  X <- matrix(0, n, 20)
  for (i in 1:10) X[i, ] <- sin(2 * pi * t) + rnorm(20, sd = 0.1)
  for (i in 11:20) X[i, ] <- 2 * sin(2 * pi * t) + rnorm(20, sd = 0.1)
  fd <- fdata(X, argvals = t)

  result <- elastic.changepoint(fd, type = "amplitude", max.iter = 3,
                                 n.mc = 50, seed = 42)

  expect_no_error(plot(result, type = "statistic"))
  expect_no_error(plot(result, type = "data"))
})
