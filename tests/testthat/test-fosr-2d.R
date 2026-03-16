# Tests for 2D function-on-scalar regression

test_that("fosr.2d returns correct structure", {
  set.seed(42)
  n <- 30
  m_s <- 8
  m_t <- 6
  m <- m_s * m_t
  argvals_s <- seq(0, 1, length.out = m_s)
  argvals_t <- seq(0, 1, length.out = m_t)

  Y <- matrix(rnorm(n * m), n, m)
  fd <- fdata(Y, argvals = seq_len(m))
  X <- cbind(rnorm(n), rnorm(n))

  fit <- fosr.2d(fd, X, argvals_s, argvals_t)

  expect_s3_class(fit, "fosr.2d")
  expect_s3_class(fit$intercept, "fdata")
  expect_s3_class(fit$beta, "fdata")
  expect_s3_class(fit$fitted, "fdata")
  expect_s3_class(fit$residuals, "fdata")
  expect_true(is.numeric(fit$r.squared))
  expect_equal(nrow(fit$fitted$data), n)
  expect_equal(fit$grid$m1, m_s)
  expect_equal(fit$grid$m2, m_t)
})

test_that("predict.fosr.2d with NULL returns fitted values", {
  set.seed(42)
  n <- 30
  m_s <- 8; m_t <- 6
  m <- m_s * m_t
  Y <- matrix(rnorm(n * m), n, m)
  fd <- fdata(Y, argvals = seq_len(m))
  X <- cbind(rnorm(n))

  fit <- fosr.2d(fd, X, seq(0, 1, length.out = m_s), seq(0, 1, length.out = m_t))

  pred_null <- predict(fit)
  expect_s3_class(pred_null, "fdata")
  expect_equal(nrow(pred_null$data), n)
})

test_that("predict.fosr.2d works with new predictors", {
  set.seed(42)
  n <- 30
  m_s <- 8; m_t <- 6
  m <- m_s * m_t
  Y <- matrix(rnorm(n * m), n, m)
  fd <- fdata(Y, argvals = seq_len(m))
  X <- cbind(rnorm(n))

  fit <- fosr.2d(fd, X, seq(0, 1, length.out = m_s), seq(0, 1, length.out = m_t))
  pred <- predict(fit, cbind(rnorm(5)))

  expect_s3_class(pred, "fdata")
  expect_equal(nrow(pred$data), 5)
})

test_that("print.fosr.2d works", {
  set.seed(42)
  n <- 30
  m_s <- 8; m_t <- 6
  m <- m_s * m_t
  Y <- matrix(rnorm(n * m), n, m)
  fd <- fdata(Y, argvals = seq_len(m))
  X <- cbind(rnorm(n))

  fit <- fosr.2d(fd, X, seq(0, 1, length.out = m_s), seq(0, 1, length.out = m_t))
  expect_output(print(fit), "2D Function-on-Scalar")
})

test_that("fosr.2d validates input", {
  expect_error(fosr.2d("not_fdata", matrix(1), 1:5, 1:5), "fdata")
})
