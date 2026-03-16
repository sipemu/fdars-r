# Extended tests for elastic regression (predict, print, plot coverage)

# Helper
.test_elastic_data <- function(n = 30, m = 20, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1, 1, 0.3)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(m, sd = 0.1)
    y[i] <- 2 * amp + rnorm(1, sd = 0.3)
  }
  fd <- fdata(X, argvals = t)
  list(fd = fd, y = y)
}

# =============================================================================
# predict.elastic.regression
# =============================================================================

test_that("predict.elastic.regression returns predictions", {
  d <- .test_elastic_data()
  fit <- elastic.regression(d$fd, d$y, max.iter = 3)

  pred <- predict(fit, d$fd[1:5, ])
  expect_equal(length(pred), 5)
  expect_true(is.numeric(pred))
  expect_true(all(is.finite(pred)))
})

test_that("predict.elastic.regression auto-wraps matrices", {
  d <- .test_elastic_data()
  fit <- elastic.regression(d$fd, d$y, max.iter = 3)

  pred <- predict(fit, d$fd$data[1:3, , drop = FALSE])
  expect_equal(length(pred), 3)
})

# =============================================================================
# predict.elastic.logistic
# =============================================================================

test_that("predict.elastic.logistic returns probabilities and classes", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 30
  X <- matrix(0, n, 20)
  for (i in 1:n) {
    amp <- if (i <= 15) rnorm(1, 1, 0.2) else rnorm(1, -1, 0.2)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(20, sd = 0.1)
  }
  fd <- fdata(X, argvals = t)
  y <- c(rep(0, 15), rep(1, 15))
  fit <- elastic.logistic(fd, y, max.iter = 3)

  probs <- predict(fit, fd[1:5, ], type = "response")
  expect_equal(length(probs), 5)
  expect_true(all(probs >= 0 & probs <= 1))

  classes <- predict(fit, fd[1:5, ], type = "class")
  expect_equal(length(classes), 5)
  expect_true(all(classes %in% c(0, 1)))
})

# =============================================================================
# print methods
# =============================================================================

test_that("print.elastic.logistic works", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 30
  X <- matrix(0, n, 20)
  for (i in 1:n) {
    amp <- if (i <= 15) 1 else -1
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(20, sd = 0.1)
  }
  fd <- fdata(X, argvals = t)
  y <- c(rep(0, 15), rep(1, 15))
  fit <- elastic.logistic(fd, y, max.iter = 3)

  expect_output(print(fit), "Logistic|Elastic|Classification", ignore.case = TRUE)
})

# =============================================================================
# plot.elastic.pcr
# =============================================================================

test_that("plot.elastic.pcr runs without error", {
  d <- .test_elastic_data()
  fit <- elastic.pcr(d$fd, d$y, ncomp = 2, max.iter = 3)

  pdf(NULL)
  expect_no_error(plot(fit, type = "fit"))
  expect_no_error(plot(fit, type = "coefficients"))
  dev.off()
})

# =============================================================================
# Classification: fclassif.cv
# =============================================================================

test_that("fclassif.cv returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 15)
  n <- 40
  X <- matrix(0, n, 15)
  y <- integer(n)
  for (i in 1:n) {
    cls <- if (i <= 20) 0L else 1L
    amp <- if (cls == 0) rnorm(1, 1, 0.2) else rnorm(1, -1, 0.2)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(15, sd = 0.1)
    y[i] <- cls
  }
  fd <- fdata(X, argvals = t)

  result <- fclassif.cv(fd, y, method = "lda", ncomp = 3, nfold = 3, seed = 42)

  expect_true(!is.null(result))
  expect_true(!is.null(result$error.rate))
})

test_that("fclassif.cv works with qda", {
  set.seed(42)
  t <- seq(0, 1, length.out = 15)
  n <- 40
  X <- matrix(0, n, 15)
  y <- integer(n)
  for (i in 1:n) {
    cls <- if (i <= 20) 0L else 1L
    amp <- if (cls == 0) rnorm(1, 1, 0.2) else rnorm(1, -1, 0.2)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(15, sd = 0.1)
    y[i] <- cls
  }
  fd <- fdata(X, argvals = t)

  result <- fclassif.cv(fd, y, method = "qda", ncomp = 3, nfold = 3, seed = 42)
  expect_true(!is.null(result))
})

test_that("fclassif.cv works with knn", {
  set.seed(42)
  t <- seq(0, 1, length.out = 15)
  n <- 40
  X <- matrix(0, n, 15)
  y <- integer(n)
  for (i in 1:n) {
    cls <- if (i <= 20) 0L else 1L
    amp <- if (cls == 0) rnorm(1, 1, 0.2) else rnorm(1, -1, 0.2)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(15, sd = 0.1)
    y[i] <- cls
  }
  fd <- fdata(X, argvals = t)

  result <- fclassif.cv(fd, y, method = "knn", ncomp = 3, nfold = 3, seed = 42)
  expect_true(!is.null(result))
})
