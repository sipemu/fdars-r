# Tests for cv.fdata (unified k-fold cross-validation)

test_that("cv.fdata regression returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 30
  X <- matrix(0, n, 20)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1, 1, 0.3)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(20, sd = 0.1)
    y[i] <- 2 * amp + rnorm(1, sd = 0.3)
  }
  fd <- fdata(X, argvals = t)

  result <- cv.fdata(fd, y,
    fit.fn = function(fd, y, ...) fregre.lm(fd, y, ncomp = 3),
    kfold = 5, seed = 42)

  expect_s3_class(result, "cv.fdata")
  expect_equal(length(result$oof.predictions), n)
  expect_true(is.numeric(result$oof.predictions))
  expect_equal(length(result$folds), n)
  expect_true(!is.null(result$metrics))
  expect_true(!is.null(result$fold.metrics))
  expect_equal(result$type, "regression")
  expect_equal(result$kfold, 5L)
})

test_that("cv.fdata classification returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 30
  X <- matrix(0, n, 20)
  for (i in 1:n) {
    amp <- if (i <= 15) rnorm(1, 1, 0.2) else rnorm(1, -1, 0.2)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(20, sd = 0.1)
  }
  fd <- fdata(X, argvals = t)
  y <- factor(c(rep("A", 15), rep("B", 15)))

  result <- cv.fdata(fd, y,
    fit.fn = function(fd, y, ...) {
      fregre.lm(fd, as.numeric(y == "B"), ncomp = 3)
    },
    predict.fn = function(model, newdata, ...) {
      p <- predict(model, newdata)
      ifelse(p > 0.5, "B", "A")
    },
    kfold = 5, seed = 42)

  expect_s3_class(result, "cv.fdata")
  expect_equal(result$type, "classification")
  expect_equal(length(result$oof.predictions), n)
})

test_that("cv.fdata with nrep > 1 returns repeated CV results", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(300), 30, 10), argvals = seq(0, 1, length.out = 10))
  y <- rnorm(30)

  result <- cv.fdata(fd, y,
    fit.fn = function(fd, y, ...) fregre.lm(fd, y, ncomp = 2),
    kfold = 3, nrep = 2, seed = 42)

  expect_s3_class(result, "cv.fdata")
  expect_true(!is.null(result$oof.matrix))
  expect_equal(ncol(result$oof.matrix), 2)
  expect_true(!is.null(result$oof.sd))
  expect_true(!is.null(result$rep.metrics))
})

test_that("cv.fdata with custom predict.fn works", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(300), 30, 10), argvals = seq(0, 1, length.out = 10))
  y <- rnorm(30)

  result <- cv.fdata(fd, y,
    fit.fn = function(fd, y, ...) fregre.lm(fd, y, ncomp = 2),
    predict.fn = function(model, newdata, ...) predict(model, newdata),
    kfold = 5, seed = 42)

  expect_equal(length(result$oof.predictions), 30)
})

test_that("cv.fdata validates input", {
  fd <- fdata(matrix(rnorm(100), 10, 10), argvals = seq(0, 1, length.out = 10))
  expect_error(cv.fdata("not_fdata", 1:10, identity), "fdata")
  expect_error(cv.fdata(fd, 1:5, identity), "Length")
  expect_error(cv.fdata(fd, 1:10, "not_a_function"), "function")
})

test_that("print.cv.fdata works", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
  y <- rnorm(20)

  result <- cv.fdata(fd, y,
    fit.fn = function(fd, y, ...) fregre.lm(fd, y, ncomp = 2),
    kfold = 3, seed = 42)

  expect_output(print(result), "Cross-Validation|cv\\.fdata|fold", ignore.case = TRUE)
})

test_that("plot.cv.fdata runs without error", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
  y <- rnorm(20)

  result <- cv.fdata(fd, y,
    fit.fn = function(fd, y, ...) fregre.lm(fd, y, ncomp = 2),
    kfold = 3, seed = 42)

  pdf(NULL)
  expect_no_error(plot(result))
  dev.off()
})
