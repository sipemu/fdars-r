test_that("fpca tolerance band returns correct structure", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
  band <- tolerance.band(fd, method = "fpca", nb = 50, seed = 42)
  expect_s3_class(band, "tolerance.band")
  expect_equal(length(band$lower), 10)
  expect_equal(length(band$upper), 10)
  expect_equal(length(band$center), 10)
  expect_equal(band$method, "fpca")
  expect_equal(band$coverage, 0.95)
})

test_that("conformal tolerance band works", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
  band <- tolerance.band(fd, method = "conformal", seed = 42)
  expect_s3_class(band, "tolerance.band")
  expect_equal(band$method, "conformal")
})

test_that("scb tolerance band works", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
  band <- tolerance.band(fd, method = "scb", nb = 50, seed = 42)
  expect_s3_class(band, "tolerance.band")
  expect_equal(band$method, "scb")
})

test_that("upper >= lower for all tolerance bands", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
  band <- tolerance.band(fd, method = "fpca", nb = 50, seed = 42)
  expect_true(all(band$upper >= band$lower))
})

test_that("tolerance band input validation works", {
  expect_error(tolerance.band("not_fdata"))
})

test_that("print method works without error", {
  set.seed(42)
  fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
  band <- tolerance.band(fd, method = "fpca", nb = 50, seed = 42)
  expect_output(print(band), "Functional Tolerance Band")
})
