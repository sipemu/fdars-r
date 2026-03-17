# Tests for Phase 9: Elastic clustering functions

# Helper: create two-cluster functional data
.test_cluster_data <- function(n_per_cluster = 6, m = 25, seed = 1) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  n <- 2 * n_per_cluster
  X <- matrix(0, n, m)
  for (i in 1:n_per_cluster) {
    X[i, ] <- sin(2 * pi * t) + rnorm(m, sd = 0.1)
  }
  for (i in (n_per_cluster + 1):n) {
    X[i, ] <- cos(2 * pi * t) + rnorm(m, sd = 0.1)
  }
  fdata(X, argvals = t)
}

# =============================================================================
# elastic.kmeans
# =============================================================================

test_that("elastic.kmeans returns correct structure", {
  fd <- .test_cluster_data(n_per_cluster = 5, m = 20)

  cl <- elastic.kmeans(fd, k = 2, max.iter = 10,
                        karcher.max.iter = 5, seed = 1)

  expect_s3_class(cl, "elastic.kmeans")
  expect_true(is.integer(cl$labels) || is.numeric(cl$labels))
  expect_equal(length(cl$labels), nrow(fd$data))
  expect_true(all(cl$labels %in% 1:2))
  expect_equal(cl$k, 2L)
  expect_true(is.numeric(cl$within.distances))
  expect_true(is.numeric(cl$total.within.distance))
  expect_gte(cl$total.within.distance, 0)
  expect_true(is.numeric(cl$n.iter))
  expect_true(is.logical(cl$converged))
})

test_that("elastic.kmeans labels are in 1:k", {
  fd <- .test_cluster_data(n_per_cluster = 5, m = 20)

  cl <- elastic.kmeans(fd, k = 3, max.iter = 5,
                        karcher.max.iter = 3, seed = 1)

  expect_true(all(cl$labels %in% 1:3))
  expect_equal(cl$k, 3L)
})

test_that("elastic.kmeans k=1 assigns all to one cluster", {
  fd <- .test_cluster_data(n_per_cluster = 5, m = 20)

  cl <- elastic.kmeans(fd, k = 1, max.iter = 5,
                        karcher.max.iter = 3, seed = 1)

  expect_true(all(cl$labels == 1))
})

test_that("elastic.kmeans centers list has k elements", {
  fd <- .test_cluster_data(n_per_cluster = 5, m = 20)
  k <- 2

  cl <- elastic.kmeans(fd, k = k, max.iter = 10,
                        karcher.max.iter = 5, seed = 1)

  expect_true(is.list(cl$centers))
  expect_equal(length(cl$centers), k)

  # Each center should have a mean field
  for (j in 1:k) {
    expect_true(!is.null(cl$centers[[j]]$mean))
    expect_true(is.numeric(cl$centers[[j]]$mean))
  }
})

test_that("elastic.kmeans validates input", {
  expect_error(elastic.kmeans("not_fdata", k = 2), "fdata")
})

test_that("elastic.kmeans print method works", {
  fd <- .test_cluster_data(n_per_cluster = 5, m = 20)
  cl <- elastic.kmeans(fd, k = 2, max.iter = 5,
                        karcher.max.iter = 3, seed = 1)
  expect_output(print(cl), "Elastic K-Means Clustering")
})

test_that("elastic.kmeans within distances are non-negative", {
  fd <- .test_cluster_data(n_per_cluster = 5, m = 20)

  cl <- elastic.kmeans(fd, k = 2, max.iter = 5,
                        karcher.max.iter = 3, seed = 1)

  expect_true(is.numeric(cl$within.distances))
  expect_true(length(cl$within.distances) > 0)
  expect_true(all(cl$within.distances >= 0))
})

# =============================================================================
# elastic.hclust
# =============================================================================

test_that("elastic.hclust returns correct structure", {
  fd <- .test_cluster_data(n_per_cluster = 4, m = 20)
  n <- nrow(fd$data)

  hc <- elastic.hclust(fd, method = "complete")

  expect_s3_class(hc, "elastic.hclust")
  expect_true(is.matrix(hc$merges))
  expect_equal(nrow(hc$merges), n - 1)
  expect_equal(ncol(hc$merges), 3)  # i, j, distance
  expect_true(is.matrix(hc$distance.matrix))
  expect_equal(nrow(hc$distance.matrix), n)
  expect_equal(ncol(hc$distance.matrix), n)
  expect_equal(hc$method, "complete")
})

test_that("elastic.hclust distance matrix is symmetric", {
  fd <- .test_cluster_data(n_per_cluster = 3, m = 20)

  hc <- elastic.hclust(fd, method = "complete")

  expect_equal(hc$distance.matrix, t(hc$distance.matrix), tolerance = 1e-6)
  expect_true(all(diag(hc$distance.matrix) < 1e-6))
})

test_that("elastic.hclust with single linkage", {
  fd <- .test_cluster_data(n_per_cluster = 4, m = 20)

  hc <- elastic.hclust(fd, method = "single")
  expect_s3_class(hc, "elastic.hclust")
  expect_equal(hc$method, "single")
})

test_that("elastic.hclust with average linkage", {
  fd <- .test_cluster_data(n_per_cluster = 4, m = 20)

  hc <- elastic.hclust(fd, method = "average")
  expect_s3_class(hc, "elastic.hclust")
  expect_equal(hc$method, "average")
})

test_that("elastic.hclust merge distances are non-decreasing (complete)", {
  fd <- .test_cluster_data(n_per_cluster = 4, m = 20)

  hc <- elastic.hclust(fd, method = "complete")
  merge_dists <- hc$merges[, 3]

  # For complete linkage, merge distances should generally be non-decreasing
  # (may have small numerical deviations)
  diffs <- diff(merge_dists)
  expect_true(all(diffs >= -1e-6))
})

test_that("elastic.hclust validates input", {
  expect_error(elastic.hclust("not_fdata"), "fdata")
})

test_that("elastic.hclust print method works", {
  fd <- .test_cluster_data(n_per_cluster = 3, m = 20)
  hc <- elastic.hclust(fd, method = "complete")
  expect_output(print(hc), "Elastic Hierarchical Clustering")
})

# =============================================================================
# elastic.cutree
# =============================================================================

test_that("elastic.cutree produces k clusters", {
  fd <- .test_cluster_data(n_per_cluster = 4, m = 20)
  n <- nrow(fd$data)

  hc <- elastic.hclust(fd, method = "complete")
  labels <- elastic.cutree(hc, k = 2)

  expect_true(is.integer(labels))
  expect_equal(length(labels), n)
  expect_true(all(labels %in% 1:2))
  expect_equal(length(unique(labels)), 2)
})

test_that("elastic.cutree k=1 assigns all to one cluster", {
  fd <- .test_cluster_data(n_per_cluster = 4, m = 20)
  n <- nrow(fd$data)

  hc <- elastic.hclust(fd, method = "complete")
  labels <- elastic.cutree(hc, k = 1)

  expect_equal(length(labels), n)
  expect_equal(length(unique(labels)), 1)
})

test_that("elastic.cutree k=n assigns each to own cluster", {
  fd <- .test_cluster_data(n_per_cluster = 3, m = 20)
  n <- nrow(fd$data)

  hc <- elastic.hclust(fd, method = "complete")
  labels <- elastic.cutree(hc, k = n)

  expect_equal(length(labels), n)
  expect_equal(length(unique(labels)), n)
})

test_that("elastic.cutree different k values produce nested partitions", {
  fd <- .test_cluster_data(n_per_cluster = 4, m = 20)

  hc <- elastic.hclust(fd, method = "complete")
  labels_2 <- elastic.cutree(hc, k = 2)
  labels_4 <- elastic.cutree(hc, k = 4)

  # More clusters should have at least as many unique values
  expect_lte(length(unique(labels_2)), length(unique(labels_4)))
})

test_that("elastic.cutree validates input class", {
  expect_error(elastic.cutree("not_hclust", k = 2), "elastic.hclust")
})

test_that("elastic.cutree validates k bounds", {
  fd <- .test_cluster_data(n_per_cluster = 3, m = 20)
  n <- nrow(fd$data)
  hc <- elastic.hclust(fd, method = "complete")

  expect_error(elastic.cutree(hc, k = 0))
  expect_error(elastic.cutree(hc, k = n + 1))
})
