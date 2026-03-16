# Extended tests for clustering: gmm.em, predict.cluster.gmm, cluster edge cases

# =============================================================================
# gmm.em (raw feature matrix EM)
# =============================================================================

test_that("gmm.em returns correct structure", {
  set.seed(42)
  # Two well-separated clusters in 3D
  features <- rbind(
    matrix(rnorm(30 * 3, mean = 0), 30, 3),
    matrix(rnorm(30 * 3, mean = 5), 30, 3)
  )
  result <- gmm.em(features, k = 2, seed = 42)

  expect_true(is.list(result))
  expect_equal(length(result$cluster), 60)
  expect_equal(length(unique(result$cluster)), 2)
  expect_true(is.matrix(result$membership))
  expect_equal(nrow(result$membership), 60)
  expect_equal(ncol(result$membership), 2)
  expect_true(is.finite(result$log.likelihood))
  expect_true(is.finite(result$bic))
  expect_true(is.finite(result$icl))
  expect_equal(length(result$weights), 2)
  expect_true(all(result$weights > 0))
  expect_equal(sum(result$weights), 1, tolerance = 1e-6)
})

test_that("gmm.em clusters separable data correctly", {
  set.seed(42)
  features <- rbind(
    matrix(rnorm(30 * 2, mean = -5), 30, 2),
    matrix(rnorm(30 * 2, mean = 5), 30, 2)
  )
  result <- gmm.em(features, k = 2, seed = 42)

  # Check that the two groups get different labels
  labels_g1 <- result$cluster[1:30]
  labels_g2 <- result$cluster[31:60]
  # Most of group 1 should share a label, most of group 2 should share a different label
  mode_g1 <- as.integer(names(sort(table(labels_g1), decreasing = TRUE))[1])
  mode_g2 <- as.integer(names(sort(table(labels_g2), decreasing = TRUE))[1])
  expect_true(mode_g1 != mode_g2)
  expect_gt(mean(labels_g1 == mode_g1), 0.8)
  expect_gt(mean(labels_g2 == mode_g2), 0.8)
})

test_that("gmm.em with diagonal covariance", {
  set.seed(42)
  features <- rbind(
    matrix(rnorm(20 * 3), 20, 3),
    matrix(rnorm(20 * 3, mean = 3), 20, 3)
  )
  result <- gmm.em(features, k = 2, cov.type = "diagonal", seed = 42)

  expect_equal(length(result$cluster), 40)
  expect_equal(length(unique(result$cluster)), 2)
})

test_that("gmm.em with spherical covariance", {
  set.seed(42)
  features <- rbind(
    matrix(rnorm(20 * 3), 20, 3),
    matrix(rnorm(20 * 3, mean = 3), 20, 3)
  )
  result <- gmm.em(features, k = 2, cov.type = "spherical", seed = 42)

  expect_equal(length(result$cluster), 40)
  expect_equal(length(unique(result$cluster)), 2)
})

test_that("gmm.em with 3 clusters", {
  set.seed(42)
  features <- rbind(
    matrix(rnorm(20 * 2, mean = 0), 20, 2),
    matrix(rnorm(20 * 2, mean = 5), 20, 2),
    matrix(rnorm(20 * 2, mean = -5), 20, 2)
  )
  result <- gmm.em(features, k = 3, seed = 42)

  expect_equal(length(result$cluster), 60)
  expect_equal(length(unique(result$cluster)), 3)
  expect_equal(ncol(result$membership), 3)
  expect_equal(length(result$weights), 3)
})

test_that("gmm.em membership rows sum to 1", {
  set.seed(42)
  features <- rbind(
    matrix(rnorm(20 * 2), 20, 2),
    matrix(rnorm(20 * 2, mean = 4), 20, 2)
  )
  result <- gmm.em(features, k = 2, seed = 42)

  row_sums <- rowSums(result$membership)
  expect_equal(row_sums, rep(1, 40), tolerance = 1e-6)
})

# =============================================================================
# predict.cluster.gmm
# =============================================================================

test_that("predict.cluster.gmm returns correct structure", {
  set.seed(42)
  m <- 20
  t <- seq(0, 1, length.out = m)
  n <- 40
  X <- rbind(
    matrix(rnorm(20 * m, mean = -1), 20, m),
    matrix(rnorm(20 * m, mean = 1), 20, m)
  )
  fd <- fdata(X, argvals = t)
  fit <- cluster.gmm(fd, k.range = 2, seed = 42)

  # Predict on new data
  X_new <- matrix(rnorm(5 * m, mean = -1), 5, m)
  fd_new <- fdata(X_new, argvals = t)
  pred <- predict(fit, fd_new)

  expect_true(is.list(pred))
  expect_equal(length(pred$cluster), 5)
  expect_true(is.matrix(pred$membership))
  expect_equal(nrow(pred$membership), 5)
})

test_that("predict.cluster.gmm validates input", {
  set.seed(42)
  m <- 20
  fd <- fdata(matrix(rnorm(40 * m), 40, m), argvals = seq(0, 1, length.out = m))
  fit <- cluster.gmm(fd, k.range = 2, seed = 42)

  expect_error(predict(fit, "not_fdata"), "fdata")
})

# =============================================================================
# cluster.gmm edge cases
# =============================================================================

test_that("cluster.gmm with different cov.type", {
  set.seed(42)
  m <- 15
  fd <- fdata(matrix(rnorm(30 * m), 30, m), argvals = seq(0, 1, length.out = m))

  fit_full <- cluster.gmm(fd, k.range = 2, cov.type = "full", seed = 42)
  fit_diag <- cluster.gmm(fd, k.range = 2, cov.type = "diagonal", seed = 42)

  expect_s3_class(fit_full, "cluster.gmm")
  expect_s3_class(fit_diag, "cluster.gmm")
})

test_that("cluster.gmm with nbasis parameter", {
  set.seed(42)
  m <- 20
  fd <- fdata(matrix(rnorm(30 * m), 30, m), argvals = seq(0, 1, length.out = m))

  fit <- cluster.gmm(fd, k.range = 2, nbasis = 5, seed = 42)
  expect_s3_class(fit, "cluster.gmm")
  expect_true(fit$d > 0)
})

# =============================================================================
# cluster.fcm edge cases
# =============================================================================

test_that("cluster.fcm with different fuzziness parameter", {
  set.seed(42)
  m <- 15
  fd <- fdata(matrix(rnorm(30 * m), 30, m), argvals = seq(0, 1, length.out = m))

  fit1 <- cluster.fcm(fd, ncl = 2, m = 1.5, seed = 42)
  fit2 <- cluster.fcm(fd, ncl = 2, m = 3.0, seed = 42)

  expect_s3_class(fit1, "cluster.fcm")
  expect_s3_class(fit2, "cluster.fcm")
  # Higher fuzziness should give more uniform membership
  max_memb_1 <- max(fit1$membership)
  max_memb_2 <- max(fit2$membership)
  expect_gt(max_memb_1, max_memb_2)
})

# =============================================================================
# cluster.init
# =============================================================================

test_that("cluster.init returns valid initial centroids", {
  set.seed(42)
  m <- 15
  fd <- fdata(matrix(rnorm(30 * m), 30, m), argvals = seq(0, 1, length.out = m))

  init <- cluster.init(fd, ncl = 3, seed = 42)

  expect_s3_class(init, "fdata")
  expect_equal(nrow(init$data), 3)
  expect_equal(ncol(init$data), m)
})
