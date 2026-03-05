test_that("detect.landmarks finds known peaks in sinusoidal data", {
  t <- seq(0, 1, length.out = 200)
  X <- matrix(sin(2 * pi * t), nrow = 1)
  fd <- fdata(X, argvals = t)
  lms <- detect.landmarks(fd, kind = "peak", min.prominence = 0.1)
  expect_type(lms, "list")
  expect_equal(length(lms), 1)
  expect_gt(nrow(lms[[1]]), 0)
  # Peak should be near t = 0.25
  expect_true(abs(lms[[1]]$position[1] - 0.25) < 0.02)
})

test_that("detect.landmarks finds valleys", {
  t <- seq(0, 1, length.out = 200)
  X <- matrix(sin(2 * pi * t), nrow = 1)
  fd <- fdata(X, argvals = t)
  lms <- detect.landmarks(fd, kind = "valley", min.prominence = 0.1)
  expect_gt(nrow(lms[[1]]), 0)
  # Valley should be near t = 0.75
  expect_true(abs(lms[[1]]$position[1] - 0.75) < 0.02)
})

test_that("detect.landmarks finds zero crossings", {
  t <- seq(0, 1, length.out = 200)
  X <- matrix(sin(2 * pi * t), nrow = 1)
  fd <- fdata(X, argvals = t)
  lms <- detect.landmarks(fd, kind = "zero")
  expect_gt(nrow(lms[[1]]), 0)
  # Zero crossing near t = 0.5
  has_half <- any(abs(lms[[1]]$position - 0.5) < 0.02)
  expect_true(has_half)
})

test_that("landmark.register aligns shifted peak data", {
  set.seed(42)
  t <- seq(0, 1, length.out = 100)
  n <- 5
  X <- matrix(0, n, 100)
  for (i in 1:n) X[i, ] <- sin(2 * pi * (t - i / 50))
  fd <- fdata(X, argvals = t)
  lr <- landmark.register(fd, kind = "peak", min.prominence = 0.5,
                          expected.count = 1)
  expect_s3_class(lr, "landmark.register")
  expect_s3_class(lr$registered, "fdata")
  expect_equal(nrow(lr$registered$data), n)
  expect_equal(ncol(lr$registered$data), 100)
  expect_type(lr$target_landmarks, "double")
})

test_that("landmark.register print method works", {
  set.seed(42)
  t <- seq(0, 1, length.out = 100)
  X <- matrix(0, 5, 100)
  for (i in 1:5) X[i, ] <- sin(2 * pi * (t - i / 50))
  fd <- fdata(X, argvals = t)
  lr <- landmark.register(fd, kind = "peak", min.prominence = 0.5,
                          expected.count = 1)
  expect_output(print(lr), "Landmark Registration")
})

test_that("detect.landmarks validates input", {
  expect_error(detect.landmarks("not_fdata"))
})
