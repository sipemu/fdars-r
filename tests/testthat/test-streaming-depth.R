test_that("streaming FM depth returns valid values", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(200), 20, 10))
  d <- streaming.depth(fd, method = "FM")
  expect_equal(length(d), 20)
  expect_true(all(is.finite(d)))
  expect_true(all(d >= 0))
})

test_that("streaming MBD depth values in [0, 1]", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(200), 20, 10))
  d <- streaming.depth(fd, method = "MBD")
  expect_equal(length(d), 20)
  expect_true(all(d >= 0 & d <= 1))
})

test_that("streaming BD depth works", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(200), 20, 10))
  d <- streaming.depth(fd, method = "BD")
  expect_equal(length(d), 20)
  expect_true(all(d >= 0 & d <= 1))
})

test_that("streaming depth with reference sample works", {
  set.seed(42)
  ref <- fdata(matrix(rnorm(200), 20, 10))
  new_data <- fdata(matrix(rnorm(50), 5, 10))
  d <- streaming.depth(new_data, ref, method = "FM")
  expect_equal(length(d), 5)
  expect_true(all(is.finite(d)))
})

test_that("depth.streaming alias works", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(200), 20, 10))
  d1 <- streaming.depth(fd, method = "FM")
  d2 <- depth.streaming(fd, method = "FM")
  expect_equal(d1, d2)
})

test_that("depth dispatcher includes streaming", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(200), 20, 10))
  d <- depth(fd, method = "streaming")
  expect_equal(length(d), 20)
  expect_true(all(is.finite(d)))
})

test_that("streaming depth input validation works", {
  expect_error(streaming.depth("not_fdata"))
  expect_error(streaming.depth(fdata(matrix(1, 1, 5)), fdataori = "not_fdata"))
})
