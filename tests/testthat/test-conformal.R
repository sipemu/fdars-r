# Tests for conformal prediction framework

# Helper: create regression test data
.test_conformal_reg_data <- function(n = 50, m = 15, seed = 42) {
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

# Helper: create classification test data
.test_conformal_class_data <- function(n = 50, m = 15, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  y <- integer(n)
  for (i in 1:n) {
    cls <- if (i <= n / 2) 0L else 1L
    amp <- if (cls == 0) rnorm(1, 1, 0.2) else rnorm(1, -1, 0.2)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(m, sd = 0.1)
    y[i] <- cls
  }
  fd <- fdata(X, argvals = t)
  list(fd = fd, y = y)
}

# =============================================================================
# Split conformal regression
# =============================================================================

test_that("conformal.fregre.lm returns valid prediction intervals", {
  d <- .test_conformal_reg_data()
  train_fd <- d$fd[1:40, ]
  test_fd <- d$fd[41:50, ]

  result <- conformal.fregre.lm(train_fd, d$y[1:40], test_fd,
                                 ncomp = 3, alpha = 0.1, seed = 42)

  expect_true(!is.null(result$predictions))
  expect_true(!is.null(result$lower))
  expect_true(!is.null(result$upper))
  expect_equal(length(result$predictions), 10)
  expect_equal(length(result$lower), 10)
  expect_equal(length(result$upper), 10)
  expect_true(all(result$upper >= result$lower))
  expect_true(result$residual_quantile > 0)
})

test_that("conformal.fregre.np returns valid prediction intervals", {
  d <- .test_conformal_reg_data()
  train_fd <- d$fd[1:40, ]
  test_fd <- d$fd[41:50, ]

  result <- conformal.fregre.np(train_fd, d$y[1:40], test_fd,
                                 alpha = 0.1, seed = 42)

  expect_true(!is.null(result$predictions))
  expect_equal(length(result$predictions), 10)
  expect_true(all(result$upper >= result$lower))
})

test_that("conformal.fregre.lm validates input", {
  expect_error(conformal.fregre.lm("not_fdata", 1:10, "not_fdata"), "fdata")
})

# =============================================================================
# Elastic conformal regression
# =============================================================================

test_that("conformal.elastic.regression returns valid intervals", {
  d <- .test_conformal_reg_data()
  train_fd <- d$fd[1:40, ]
  test_fd <- d$fd[41:50, ]

  result <- conformal.elastic.regression(train_fd, d$y[1:40], test_fd,
                                          ncomp.beta = 5, alpha = 0.1, seed = 42)

  expect_true(!is.null(result$predictions))
  expect_equal(length(result$predictions), 10)
  expect_true(all(result$upper >= result$lower))
})

test_that("conformal.elastic.pcr returns valid intervals", {
  d <- .test_conformal_reg_data()
  train_fd <- d$fd[1:40, ]
  test_fd <- d$fd[41:50, ]

  result <- conformal.elastic.pcr(train_fd, d$y[1:40], test_fd,
                                   ncomp = 3, alpha = 0.1, seed = 42)

  expect_true(!is.null(result$predictions))
  expect_equal(length(result$predictions), 10)
  expect_true(all(result$upper >= result$lower))
})

# =============================================================================
# Conformal classification
# =============================================================================

test_that("conformal.classif returns prediction sets", {
  d <- .test_conformal_class_data()
  train_fd <- d$fd[1:40, ]
  test_fd <- d$fd[41:50, ]

  result <- conformal.classif(train_fd, d$y[1:40], test_fd,
                               ncomp = 3, classifier = "lda",
                               alpha = 0.1, seed = 42)

  expect_true(!is.null(result$predicted_classes))
  expect_true(!is.null(result$set_sizes))
  expect_equal(length(result$predicted_classes), 10)
  expect_equal(length(result$set_sizes), 10)
  expect_true(result$average_set_size > 0)
})

test_that("conformal.classif works with knn classifier", {
  d <- .test_conformal_class_data()
  train_fd <- d$fd[1:40, ]
  test_fd <- d$fd[41:50, ]

  result <- conformal.classif(train_fd, d$y[1:40], test_fd,
                               ncomp = 3, classifier = "knn", k.nn = 3,
                               alpha = 0.1, seed = 42)

  expect_true(!is.null(result$predicted_classes))
  expect_equal(length(result$predicted_classes), 10)
})

test_that("conformal.logistic returns prediction sets", {
  d <- .test_conformal_class_data()
  train_fd <- d$fd[1:40, ]
  test_fd <- d$fd[41:50, ]

  result <- conformal.logistic(train_fd, d$y[1:40], test_fd,
                                ncomp = 3, alpha = 0.1, seed = 42)

  expect_true(!is.null(result$predicted_classes))
  expect_equal(length(result$predicted_classes), 10)
})

test_that("conformal.elastic.logistic returns prediction sets", {
  d <- .test_conformal_class_data()
  train_fd <- d$fd[1:40, ]
  test_fd <- d$fd[41:50, ]

  result <- conformal.elastic.logistic(train_fd, d$y[1:40], test_fd,
                                        alpha = 0.1, seed = 42)

  expect_true(!is.null(result$predicted_classes))
  expect_equal(length(result$predicted_classes), 10)
})

# =============================================================================
# CV-conformal and jackknife+
# =============================================================================

test_that("cv.conformal.regression returns valid intervals", {
  d <- .test_conformal_reg_data()
  train_fd <- d$fd[1:40, ]
  test_fd <- d$fd[41:50, ]

  result <- cv.conformal.regression(train_fd, d$y[1:40], test_fd,
                                     method = "fregre.lm", ncomp = 3,
                                     n.folds = 3, alpha = 0.1, seed = 42)

  expect_true(!is.null(result$predictions))
  expect_equal(length(result$predictions), 10)
  expect_true(all(result$upper >= result$lower))
})

test_that("cv.conformal.classification returns prediction sets", {
  d <- .test_conformal_class_data()
  train_fd <- d$fd[1:40, ]
  test_fd <- d$fd[41:50, ]

  result <- cv.conformal.classification(train_fd, d$y[1:40], test_fd,
                                          ncomp = 3, classifier = "lda",
                                          n.folds = 3, alpha = 0.1, seed = 42)

  expect_true(!is.null(result$predicted_classes))
  expect_equal(length(result$predicted_classes), 10)
})

test_that("jackknife.plus returns valid intervals", {
  d <- .test_conformal_reg_data(n = 30)
  train_fd <- d$fd[1:25, ]
  test_fd <- d$fd[26:30, ]

  result <- jackknife.plus(train_fd, d$y[1:25], test_fd,
                            method = "fregre.lm", ncomp = 3, alpha = 0.1)

  expect_true(!is.null(result$predictions))
  expect_equal(length(result$predictions), 5)
  expect_true(all(result$upper >= result$lower))
})

# =============================================================================
# Generic conformal
# =============================================================================

test_that("conformal.generic.regression works with pre-fitted model", {
  d <- .test_conformal_reg_data()
  model <- fregre.lm(d$fd, d$y, ncomp = 3)
  test_fd <- d$fd[1:10, ]

  result <- conformal.generic.regression(model, d$fd, d$y, test_fd,
                                          alpha = 0.1, seed = 42)

  expect_true(!is.null(result$predictions))
  expect_equal(length(result$predictions), 10)
  expect_true(all(result$upper >= result$lower))
})

test_that("conformal.generic.classification works with pre-fitted model", {
  d <- .test_conformal_class_data()
  model <- functional.logistic(d$fd, d$y, ncomp = 3)
  test_fd <- d$fd[1:10, ]

  result <- conformal.generic.classification(model, d$fd, d$y, test_fd,
                                              alpha = 0.1, seed = 42)

  expect_true(!is.null(result$predicted_classes))
  expect_equal(length(result$predicted_classes), 10)
})

# =============================================================================
# Wider alpha gives wider intervals
# =============================================================================

test_that("smaller alpha gives wider intervals", {
  d <- .test_conformal_reg_data()
  train_fd <- d$fd[1:40, ]
  test_fd <- d$fd[41:50, ]

  r_wide <- conformal.fregre.lm(train_fd, d$y[1:40], test_fd,
                                  ncomp = 3, alpha = 0.05, seed = 42)
  r_narrow <- conformal.fregre.lm(train_fd, d$y[1:40], test_fd,
                                    ncomp = 3, alpha = 0.2, seed = 42)

  width_wide <- mean(r_wide$upper - r_wide$lower)
  width_narrow <- mean(r_narrow$upper - r_narrow$lower)
  expect_true(width_wide >= width_narrow)
})
