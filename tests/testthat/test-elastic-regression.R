# Tests for elastic regression, logistic, and PCR

test_that("elastic.regression returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 20
  X <- matrix(0, n, 20)
  y <- numeric(n)
  for (i in 1:n) {
    phase <- runif(1, -0.1, 0.1)
    amp <- rnorm(1, 1, 0.3)
    X[i, ] <- amp * sin(2 * pi * (t + phase))
    y[i] <- 2 * amp + rnorm(1, sd = 0.1)
  }
  fd <- fdata(X, argvals = t)

  result <- elastic.regression(fd, y, ncomp.beta = 5, max.iter = 5)

  expect_s3_class(result, "elastic.regression")
  expect_true(is.finite(result$alpha))
  expect_s3_class(result$beta, "fdata")
  expect_equal(length(result$fitted.values), n)
  expect_equal(length(result$residuals), n)
  expect_true(is.finite(result$r.squared))
  expect_gte(result$r.squared, 0)
  expect_true(result$n.iter >= 1)
})

test_that("elastic.logistic returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 30
  X <- matrix(0, n, 20)
  y <- rep(c(0, 1), each = n / 2)
  for (i in 1:n) {
    if (y[i] == 0) {
      X[i, ] <- sin(2 * pi * t) + rnorm(20, sd = 0.1)
    } else {
      X[i, ] <- cos(2 * pi * t) + rnorm(20, sd = 0.1)
    }
  }
  fd <- fdata(X, argvals = t)

  result <- elastic.logistic(fd, y, ncomp.beta = 5, max.iter = 5)

  expect_s3_class(result, "elastic.logistic")
  expect_true(is.finite(result$alpha))
  expect_s3_class(result$beta, "fdata")
  expect_equal(length(result$probabilities), n)
  expect_equal(length(result$predicted.classes), n)
  expect_true(all(result$probabilities >= 0 & result$probabilities <= 1))
  expect_true(is.finite(result$accuracy))
  expect_gte(result$accuracy, 0)
  expect_lte(result$accuracy, 1)
})

test_that("elastic.pcr returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 20
  X <- matrix(0, n, 20)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1, 1, 0.3)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(20, sd = 0.1)
    y[i] <- 2 * amp + rnorm(1, sd = 0.2)
  }
  fd <- fdata(X, argvals = t)

  result <- elastic.pcr(fd, y, ncomp = 2, pca.method = "vertical", max.iter = 5)

  expect_s3_class(result, "elastic.pcr")
  expect_true(is.finite(result$alpha))
  expect_equal(length(result$coefficients), 2)
  expect_equal(length(result$fitted.values), n)
  expect_true(is.finite(result$r.squared))
  expect_gte(result$r.squared, 0)
  expect_equal(result$pca.method, "vertical")
  expect_equal(result$ncomp, 2)
})

test_that("elastic regression input validation works", {
  expect_error(elastic.regression("not_fdata", 1:10))
  expect_error(elastic.logistic("not_fdata", c(0, 1)))
  expect_error(elastic.pcr("not_fdata", 1:10))

  fd <- fdata(matrix(rnorm(50), 5, 10), argvals = seq(0, 1, length.out = 10))
  expect_error(elastic.regression(fd, 1:3))  # length mismatch
  expect_error(elastic.logistic(fd, c(0, 1)))  # length mismatch
  expect_error(elastic.pcr(fd, 1:3))  # length mismatch
})

test_that("print methods work for elastic regression", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 15
  X <- matrix(0, n, 20)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1, 1, 0.3)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(20, sd = 0.1)
    y[i] <- 2 * amp + rnorm(1, sd = 0.2)
  }
  fd <- fdata(X, argvals = t)

  result <- elastic.regression(fd, y, ncomp.beta = 5, max.iter = 5)
  expect_output(print(result), "Elastic Scalar-on-Function Regression")

  result_pcr <- elastic.pcr(fd, y, ncomp = 2, max.iter = 5)
  expect_output(print(result_pcr), "Elastic Principal Component Regression")
})

test_that("plot methods work for elastic regression", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 15
  X <- matrix(0, n, 20)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1, 1, 0.3)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(20, sd = 0.1)
    y[i] <- 2 * amp + rnorm(1, sd = 0.2)
  }
  fd <- fdata(X, argvals = t)

  result <- elastic.regression(fd, y, ncomp.beta = 5, max.iter = 5)
  expect_no_error(plot(result, type = "fit"))
  expect_no_error(plot(result, type = "beta"))
  expect_no_error(plot(result, type = "residuals"))
})
