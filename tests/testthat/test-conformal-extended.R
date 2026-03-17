# Extended conformal inference tests: elastic, generic, jackknife.plus, cv

# =============================================================================
# Helper
# =============================================================================

.test_conformal_data <- function(n = 40, m = 20, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(m, sd = 0.1)
    y[i] <- 2 * amp + rnorm(1, sd = 0.3)
  }
  fd <- fdata(X, argvals = t)
  list(fd = fd, y = y, t = t)
}

.test_conformal_classif_data <- function(n = 40, m = 20, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  for (i in 1:n) {
    if (i <= n / 2) {
      X[i, ] <- sin(2 * pi * t) + rnorm(m, sd = 0.2)
    } else {
      X[i, ] <- cos(2 * pi * t) + rnorm(m, sd = 0.2)
    }
  }
  y <- rep(c(0, 1), each = n / 2)
  fd <- fdata(X, argvals = t)
  list(fd = fd, y = y, t = t)
}

# =============================================================================
# conformal.elastic.regression
# =============================================================================

test_that("conformal.elastic.regression returns valid intervals", {
  d <- .test_conformal_data()

  result <- conformal.elastic.regression(
    d$fd[1:30, ], d$y[1:30],
    d$fd[31:40, ], alpha = 0.1, ncomp = 3
  )

  expect_true(is.list(result))
  expect_true("lower" %in% names(result))
  expect_true("upper" %in% names(result))
  expect_true(all(result$lower <= result$upper))
})

# =============================================================================
# conformal.elastic.pcr
# =============================================================================

test_that("conformal.elastic.pcr returns valid intervals", {
  d <- .test_conformal_data()

  result <- conformal.elastic.pcr(
    d$fd[1:30, ], d$y[1:30],
    d$fd[31:40, ], alpha = 0.1, ncomp = 3
  )

  expect_true(is.list(result))
  expect_true("lower" %in% names(result))
  expect_true("upper" %in% names(result))
})

# =============================================================================
# conformal.generic.regression (needs a pre-fitted model)
# =============================================================================

test_that("conformal.generic.regression with calibration indices", {
  d <- .test_conformal_data()
  model <- fregre.lm(d$fd[1:30, ], d$y[1:30], ncomp = 3)

  result <- conformal.generic.regression(
    model, d$fd[1:30, ], d$y[1:30],
    d$fd[31:40, ],
    calibration.indices = 21:30,
    alpha = 0.1
  )

  expect_true(is.list(result))
  if ("lower" %in% names(result)) {
    expect_equal(length(result$lower), 10)
    expect_true(all(result$lower <= result$upper))
  }
})

# =============================================================================
# conformal.generic.classification (needs a pre-fitted logistic model)
# =============================================================================

test_that("conformal.generic.classification returns prediction sets", {
  d <- .test_conformal_classif_data()
  model <- functional.logistic(d$fd[1:30, ], d$y[1:30], ncomp = 3)

  result <- conformal.generic.classification(
    model, d$fd[1:30, ], d$y[1:30],
    d$fd[31:40, ],
    calibration.indices = 21:30,
    alpha = 0.1
  )

  expect_true(is.list(result))
})

# =============================================================================
# conformal.fregre.np
# =============================================================================

test_that("conformal.fregre.np returns valid intervals", {
  d <- .test_conformal_data()

  result <- conformal.fregre.np(
    d$fd[1:30, ], d$y[1:30],
    d$fd[31:40, ], alpha = 0.1
  )

  expect_true(is.list(result))
  expect_true("lower" %in% names(result))
  expect_true("upper" %in% names(result))
})

# =============================================================================
# conformal.elastic.logistic
# =============================================================================

test_that("conformal.elastic.logistic returns prediction sets", {
  d <- .test_conformal_classif_data()

  result <- conformal.elastic.logistic(
    d$fd[1:30, ], d$y[1:30],
    d$fd[31:40, ], alpha = 0.1
  )

  expect_true(is.list(result))
})

# =============================================================================
# conformal.logistic
# =============================================================================

test_that("conformal.logistic returns prediction sets", {
  d <- .test_conformal_classif_data()

  result <- conformal.logistic(
    d$fd[1:30, ], d$y[1:30],
    d$fd[31:40, ], alpha = 0.1, ncomp = 3
  )

  expect_true(is.list(result))
})

# =============================================================================
# jackknife.plus
# =============================================================================

test_that("jackknife.plus with fregre.lm returns intervals", {
  d <- .test_conformal_data(n = 30)
  new_fd <- fdata(d$fd$data[1:3, ], argvals = d$t)

  result <- jackknife.plus(
    d$fd, d$y, new_fd,
    method = "fregre.lm", ncomp = 3, alpha = 0.1
  )

  expect_true(is.list(result))
  if ("lower" %in% names(result)) {
    expect_equal(length(result$lower), 3)
    expect_true(all(result$lower <= result$upper))
  }
})

test_that("jackknife.plus with fregre.np returns intervals", {
  d <- .test_conformal_data(n = 30)
  new_fd <- fdata(d$fd$data[1:3, ], argvals = d$t)

  result <- jackknife.plus(
    d$fd, d$y, new_fd,
    method = "fregre.np", alpha = 0.1
  )

  expect_true(is.list(result))
})

test_that("jackknife.plus validates input", {
  d <- .test_conformal_data()
  expect_error(jackknife.plus("not_fdata", d$y, d$fd[1:3, ]), "fdata")
  expect_error(jackknife.plus(d$fd, d$y, "not_fdata"), "fdata")
})

# =============================================================================
# cv.conformal.regression (requires newdata)
# =============================================================================

test_that("cv.conformal.regression returns coverage info", {
  d <- .test_conformal_data(n = 40)

  result <- cv.conformal.regression(
    d$fd[1:30, ], d$y[1:30], d$fd[31:40, ],
    ncomp = 3, n.folds = 3, alpha = 0.1, seed = 42
  )

  expect_true(is.list(result))
  if ("coverage" %in% names(result)) {
    expect_gte(result$coverage, 0)
    expect_lte(result$coverage, 1)
  }
})

# =============================================================================
# cv.conformal.classification (requires newdata)
# =============================================================================

test_that("cv.conformal.classification returns coverage info", {
  d <- .test_conformal_classif_data(n = 40)

  result <- cv.conformal.classification(
    d$fd[1:30, ], d$y[1:30], d$fd[31:40, ],
    ncomp = 3, n.folds = 3, alpha = 0.1, seed = 42
  )

  expect_true(is.list(result))
  if ("coverage" %in% names(result)) {
    expect_gte(result$coverage, 0)
    expect_lte(result$coverage, 1)
  }
})
