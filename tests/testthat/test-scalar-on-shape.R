# Tests for scalar-on-shape regression

# Helper: create test data
.test_sos_data <- function(n = 40, m = 20, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1, 1, 0.3)
    phase <- rnorm(1, 0, 0.1)
    X[i, ] <- amp * sin(2 * pi * (t + phase)) + rnorm(m, sd = 0.05)
    y[i] <- 2 * amp + rnorm(1, sd = 0.3)
  }
  fd <- fdata(X, argvals = t)
  list(fd = fd, y = y)
}

# =============================================================================
# scalar.on.shape
# =============================================================================

test_that("scalar.on.shape returns correct structure", {
  d <- .test_sos_data()
  fit <- scalar.on.shape(d$fd, d$y, nbasis = 5, max.iter.outer = 3)

  expect_s3_class(fit, "scalar.on.shape")
  expect_equal(length(fit$fitted.values), 40)
  expect_equal(length(fit$residuals), 40)
  expect_true(is.numeric(fit$beta))
  expect_equal(length(fit$beta), 20)
  expect_true(is.numeric(fit$r.squared))
  expect_true(is.numeric(fit$r.squared) && is.finite(fit$r.squared))
  expect_true(is.numeric(fit$sse))
  expect_true(is.numeric(fit$shape.scores))
  expect_equal(length(fit$shape.scores), 40)
  expect_true(!is.null(fit$gammas))
  expect_true(!is.null(fit$h.coefficients))
  expect_true(!is.null(fit$g.coefficients))
  expect_equal(fit$index.method, "identity")
})

test_that("scalar.on.shape works with polynomial index", {
  d <- .test_sos_data()
  fit <- scalar.on.shape(d$fd, d$y, nbasis = 5,
                          index.method = "polynomial", index.param = 3,
                          max.iter.outer = 3)
  expect_s3_class(fit, "scalar.on.shape")
  expect_equal(fit$index.method, "polynomial")
})

test_that("scalar.on.shape works with nadaraya.watson index", {
  d <- .test_sos_data()
  fit <- scalar.on.shape(d$fd, d$y, nbasis = 5,
                          index.method = "nadaraya.watson",
                          max.iter.outer = 3)
  expect_s3_class(fit, "scalar.on.shape")
  expect_equal(fit$index.method, "nadaraya.watson")
})

test_that("scalar.on.shape validates input", {
  expect_error(scalar.on.shape("not_fdata", 1:10), "fdata")
})

# =============================================================================
# predict.scalar.on.shape
# =============================================================================

test_that("predict.scalar.on.shape returns predictions", {
  d <- .test_sos_data()
  fit <- scalar.on.shape(d$fd, d$y, nbasis = 5, max.iter.outer = 3)

  pred <- predict(fit, d$fd[1:5, ])
  expect_equal(length(pred), 5)
  expect_true(is.numeric(pred))
  expect_true(all(is.finite(pred)))
})

test_that("predict.scalar.on.shape validates input", {
  d <- .test_sos_data()
  fit <- scalar.on.shape(d$fd, d$y, nbasis = 5, max.iter.outer = 3)
  expect_error(predict(fit, "not_fdata"), "fdata")
})

# =============================================================================
# print.scalar.on.shape
# =============================================================================

test_that("print.scalar.on.shape works", {
  d <- .test_sos_data()
  fit <- scalar.on.shape(d$fd, d$y, nbasis = 5, max.iter.outer = 3)
  expect_output(print(fit), "Scalar-on-Shape")
})
