# Tests for fdata arithmetic operations and edge cases

# =============================================================================
# Ops.fdata: arithmetic operators
# =============================================================================

test_that("fdata + fdata works correctly", {
  t <- seq(0, 1, length.out = 10)
  fd1 <- fdata(matrix(1:10, 1, 10), argvals = t)
  fd2 <- fdata(matrix(10:1, 1, 10), argvals = t)

  result <- fd1 + fd2
  expect_s3_class(result, "fdata")
  expect_equal(as.vector(result$data), rep(11, 10))
  expect_equal(result$argvals, t)
})

test_that("fdata - fdata works correctly", {
  t <- seq(0, 1, length.out = 10)
  fd1 <- fdata(matrix(rep(5, 10), 1, 10), argvals = t)
  fd2 <- fdata(matrix(rep(3, 10), 1, 10), argvals = t)

  result <- fd1 - fd2
  expect_equal(as.vector(result$data), rep(2, 10))
})

test_that("fdata * fdata works correctly", {
  t <- seq(0, 1, length.out = 10)
  fd1 <- fdata(matrix(rep(3, 10), 1, 10), argvals = t)
  fd2 <- fdata(matrix(rep(4, 10), 1, 10), argvals = t)

  result <- fd1 * fd2
  expect_equal(as.vector(result$data), rep(12, 10))
})

test_that("fdata / fdata works correctly", {
  t <- seq(0, 1, length.out = 10)
  fd1 <- fdata(matrix(rep(12, 10), 1, 10), argvals = t)
  fd2 <- fdata(matrix(rep(4, 10), 1, 10), argvals = t)

  result <- fd1 / fd2
  expect_equal(as.vector(result$data), rep(3, 10))
})

test_that("fdata ^ scalar works correctly", {
  t <- seq(0, 1, length.out = 10)
  fd <- fdata(matrix(rep(3, 10), 1, 10), argvals = t)

  result <- fd^2
  expect_equal(as.vector(result$data), rep(9, 10))
})

test_that("fdata + scalar works correctly", {
  t <- seq(0, 1, length.out = 10)
  fd <- fdata(matrix(1:10, 1, 10), argvals = t)

  result <- fd + 100
  expect_equal(as.vector(result$data), 101:110)
})

test_that("scalar + fdata works correctly", {
  t <- seq(0, 1, length.out = 10)
  fd <- fdata(matrix(1:10, 1, 10), argvals = t)

  result <- 100 + fd
  expect_equal(as.vector(result$data), 101:110)
})

test_that("fdata * scalar works correctly", {
  t <- seq(0, 1, length.out = 10)
  fd <- fdata(matrix(rep(5, 10), 1, 10), argvals = t)

  result <- fd * 3
  expect_equal(as.vector(result$data), rep(15, 10))
})

test_that("fdata / scalar works correctly", {
  t <- seq(0, 1, length.out = 10)
  fd <- fdata(matrix(rep(12, 10), 1, 10), argvals = t)

  result <- fd / 4
  expect_equal(as.vector(result$data), rep(3, 10))
})

test_that("fdata - scalar works correctly", {
  t <- seq(0, 1, length.out = 10)
  fd <- fdata(matrix(rep(10, 10), 1, 10), argvals = t)

  result <- fd - 3
  expect_equal(as.vector(result$data), rep(7, 10))
})

test_that("Ops.fdata rejects mismatched argvals", {
  fd1 <- fdata(matrix(1:10, 1, 10), argvals = seq(0, 1, length.out = 10))
  fd2 <- fdata(matrix(1:10, 1, 10), argvals = seq(0, 2, length.out = 10))

  expect_error(fd1 + fd2, "argvals")
})

test_that("Ops.fdata rejects mismatched dimensions", {
  t <- seq(0, 1, length.out = 10)
  fd1 <- fdata(matrix(1:10, 1, 10), argvals = t)
  fd2 <- fdata(matrix(1:20, 2, 10), argvals = t)

  expect_error(fd1 + fd2, "dimensions")
})

test_that("Ops.fdata preserves argvals and names", {
  t <- seq(0, 1, length.out = 10)
  fd <- fdata(matrix(1:10, 1, 10), argvals = t,
              names = list(main = "test", xlab = "Time", ylab = "Value"))

  result <- fd * 2
  expect_equal(result$argvals, t)
  expect_equal(result$names$main, "test")
})

test_that("Ops.fdata works with multi-row fdata", {
  t <- seq(0, 1, length.out = 10)
  fd1 <- fdata(matrix(rnorm(30), 3, 10), argvals = t)
  fd2 <- fdata(matrix(rnorm(30), 3, 10), argvals = t)

  result <- fd1 + fd2
  expect_s3_class(result, "fdata")
  expect_equal(nrow(result$data), 3)
  expect_equal(ncol(result$data), 10)
  expect_equal(result$data, fd1$data + fd2$data)
})

# =============================================================================
# Comparison operators
# =============================================================================

test_that("fdata > scalar works", {
  t <- seq(0, 1, length.out = 5)
  fd <- fdata(matrix(c(1, 3, 5, 2, 4), 1, 5), argvals = t)

  result <- fd > 3
  expect_s3_class(result, "fdata")
  expect_equal(as.vector(result$data), c(FALSE, FALSE, TRUE, FALSE, TRUE))
})

test_that("fdata == fdata works", {
  t <- seq(0, 1, length.out = 5)
  fd1 <- fdata(matrix(c(1, 2, 3, 4, 5), 1, 5), argvals = t)
  fd2 <- fdata(matrix(c(1, 0, 3, 0, 5), 1, 5), argvals = t)

  result <- fd1 == fd2
  expect_equal(as.vector(result$data), c(TRUE, FALSE, TRUE, FALSE, TRUE))
})

# =============================================================================
# fdata subsetting edge cases
# =============================================================================

test_that("[.fdata single row returns fdata", {
  t <- seq(0, 1, length.out = 10)
  fd <- fdata(matrix(rnorm(30), 3, 10), argvals = t)

  sub <- fd[1, ]
  expect_s3_class(sub, "fdata")
  expect_equal(nrow(sub$data), 1)
  expect_equal(ncol(sub$data), 10)
})

test_that("[.fdata range returns fdata", {
  t <- seq(0, 1, length.out = 10)
  fd <- fdata(matrix(rnorm(50), 5, 10), argvals = t)

  sub <- fd[2:4, ]
  expect_s3_class(sub, "fdata")
  expect_equal(nrow(sub$data), 3)
})

test_that("[.fdata with logical index", {
  t <- seq(0, 1, length.out = 10)
  fd <- fdata(matrix(rnorm(40), 4, 10), argvals = t)

  idx <- c(TRUE, FALSE, TRUE, FALSE)
  sub <- fd[idx, ]
  expect_s3_class(sub, "fdata")
  expect_equal(nrow(sub$data), 2)
})
