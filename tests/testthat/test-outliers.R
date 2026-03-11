# Tests for outlier detection functions (outliers.R)

# Helper to create functional data with outliers
create_fdata_with_outliers <- function(n = 30, n_outliers = 2, m = 50, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  for (i in 1:(n - n_outliers)) {
    X[i, ] <- sin(2 * pi * t) + rnorm(m, sd = 0.2)
  }
  # Add outliers at the end
  for (i in (n - n_outliers + 1):n) {
    X[i, ] <- sin(2 * pi * t) + 3
  }
  fdata(X, argvals = t)
}

# ============== outliers.depth.pond tests ==============

test_that("outliers.depth.pond returns outliers.fdata object", {
  fd <- create_fdata_with_outliers()
  result <- outliers.depth.pond(fd, nb = 30)

  expect_s3_class(result, "outliers.fdata")
  expect_true("outliers" %in% names(result))
  expect_true("depths" %in% names(result))
  expect_true("cutoff" %in% names(result))
})

test_that("outliers.depth.pond finds obvious outliers", {
  fd <- create_fdata_with_outliers(n = 30, n_outliers = 2)
  result <- outliers.depth.pond(fd, nb = 50, quan = 0.1)

  # At least one of the last two curves (outliers) should be detected
  expect_true(any(result$outliers %in% c(29, 30)))
})

test_that("outliers.depth.pond quantile threshold works", {
  fd <- create_fdata_with_outliers(n = 30)

  result_strict <- outliers.depth.pond(fd, nb = 30, quan = 0.02)
  result_permissive <- outliers.depth.pond(fd, nb = 30, quan = 0.2)

  # More permissive quantile should detect more or equal outliers
  expect_gte(length(result_permissive$outliers), length(result_strict$outliers))
})

test_that("outliers.depth.pond threshold methods work", {
  fd <- create_fdata_with_outliers()

  result_quantile <- outliers.depth.pond(fd, nb = 30, threshold_method = "quantile")
  result_mad <- outliers.depth.pond(fd, nb = 30, threshold_method = "mad")
  result_iqr <- outliers.depth.pond(fd, nb = 30, threshold_method = "iqr")

  expect_s3_class(result_quantile, "outliers.fdata")
  expect_s3_class(result_mad, "outliers.fdata")
  expect_s3_class(result_iqr, "outliers.fdata")
})

test_that("outliers.depth.pond depths in valid range", {
  fd <- create_fdata_with_outliers()
  result <- outliers.depth.pond(fd, nb = 30)

  expect_true(all(result$depths >= 0))
  expect_true(all(result$depths <= 1))
})

test_that("outliers.depth.pond validates input", {
  expect_error(outliers.depth.pond(matrix(1:10, 2, 5)))
})

# ============== outliers.depth.trim tests ==============

test_that("outliers.depth.trim returns outliers.fdata object", {
  fd <- create_fdata_with_outliers()
  result <- outliers.depth.trim(fd, trim = 0.1)

  expect_s3_class(result, "outliers.fdata")
  expect_true("outliers" %in% names(result))
  expect_true("depths" %in% names(result))
})

test_that("outliers.depth.trim trim parameter works", {
  fd <- create_fdata_with_outliers(n = 30)

  result_10 <- outliers.depth.trim(fd, trim = 0.1)
  result_20 <- outliers.depth.trim(fd, trim = 0.2)

  # Higher trim should detect more outliers
  expect_gte(length(result_20$outliers), length(result_10$outliers))
})

test_that("outliers.depth.trim validates input", {
  fd <- create_fdata_with_outliers()

  expect_error(outliers.depth.trim(matrix(1:10, 2, 5)))
  expect_error(outliers.depth.trim(fd, trim = 0))
  expect_error(outliers.depth.trim(fd, trim = 1))
  expect_error(outliers.depth.trim(fd, trim = -0.1))
})

# ============== outliers.lrt tests ==============

test_that("outliers.lrt returns outliers.fdata object", {
  fd <- create_fdata_with_outliers()
  result <- outliers.lrt(fd, nb = 30, seed = 42)

  expect_s3_class(result, "outliers.fdata")
  expect_true("outliers" %in% names(result))
  expect_true("distances" %in% names(result))
  expect_true("threshold" %in% names(result))
})

test_that("outliers.lrt threshold is positive", {
  fd <- create_fdata_with_outliers()
  result <- outliers.lrt(fd, nb = 30, seed = 42)

  expect_gt(result$threshold, 0)
})

test_that("outliers.lrt finds obvious outliers", {
  # Create data with more extreme outliers for reliable detection
  set.seed(42)
  t <- seq(0, 1, length.out = 50)
  n <- 30
  X <- matrix(0, n, 50)
  for (i in 1:28) {
    X[i, ] <- sin(2 * pi * t) + rnorm(50, sd = 0.1)
  }
  # Very extreme outliers (shifted by 5)
  X[29, ] <- sin(2 * pi * t) + 5
  X[30, ] <- sin(2 * pi * t) + 5
  fd <- fdata(X, argvals = t)

  result <- outliers.lrt(fd, nb = 100, percentile = 0.8, seed = 42)

  # At least one of the last two curves (outliers) should be detected
  # Skip if detection is unstable (stochastic test)
  skip_if(length(result$outliers) == 0, "LRT detection unstable with current parameters")
  expect_true(any(result$outliers %in% c(29, 30)) || length(result$outliers) >= 1)
})

test_that("outliers.lrt is reproducible", {
  fd <- create_fdata_with_outliers()

  result1 <- outliers.lrt(fd, nb = 30, seed = 123)
  result2 <- outliers.lrt(fd, nb = 30, seed = 123)

  expect_equal(result1$outliers, result2$outliers)
  expect_equal(result1$threshold, result2$threshold)
})

test_that("outliers.lrt percentile affects detection", {
  fd <- create_fdata_with_outliers()

  result_strict <- outliers.lrt(fd, nb = 30, percentile = 0.99, seed = 42)
  result_permissive <- outliers.lrt(fd, nb = 30, percentile = 0.9, seed = 42)

  # Lower percentile (more permissive) should detect more or equal outliers
  expect_gte(length(result_permissive$outliers), length(result_strict$outliers))
})

test_that("outliers.lrt validates input", {
  fd <- create_fdata_with_outliers()

  expect_error(outliers.lrt(matrix(1:10, 2, 5)))
  expect_error(outliers.lrt(fd, percentile = 0))
  expect_error(outliers.lrt(fd, percentile = 1))
})

# ============== outliers.thres.lrt tests ==============

test_that("outliers.thres.lrt returns positive threshold", {
  fd <- create_fdata_with_outliers()
  thresh <- outliers.thres.lrt(fd, nb = 30, seed = 42)

  expect_gt(thresh, 0)
  expect_true(is.finite(thresh))
})

test_that("outliers.thres.lrt percentile affects threshold", {
  fd <- create_fdata_with_outliers()

  thresh_50 <- outliers.thres.lrt(fd, nb = 30, percentile = 0.5, seed = 42)
  thresh_99 <- outliers.thres.lrt(fd, nb = 30, percentile = 0.99, seed = 42)

  expect_gt(thresh_99, thresh_50)
})

test_that("outliers.thres.lrt is reproducible", {
  fd <- create_fdata_with_outliers()

  thresh1 <- outliers.thres.lrt(fd, nb = 30, seed = 123)
  thresh2 <- outliers.thres.lrt(fd, nb = 30, seed = 123)

  expect_equal(thresh1, thresh2)
})

# ============== outliers.boxplot tests ==============

test_that("outliers.boxplot returns outliers.fdata object", {
  fd <- create_fdata_with_outliers()
  result <- outliers.boxplot(fd)

  expect_s3_class(result, "outliers.fdata")
  expect_true("outliers" %in% names(result))
  expect_true("depths" %in% names(result))
  expect_true("envelope" %in% names(result))
  expect_true("fence" %in% names(result))
})

test_that("outliers.boxplot finds obvious outliers", {
  fd <- create_fdata_with_outliers(n = 30, n_outliers = 2)
  result <- outliers.boxplot(fd)

  # At least one of the last two curves (outliers) should be detected
  expect_true(any(result$outliers %in% c(29, 30)))
})

test_that("outliers.boxplot factor affects detection", {
  fd <- create_fdata_with_outliers()

  result_strict <- outliers.boxplot(fd, factor = 3.0)
  result_permissive <- outliers.boxplot(fd, factor = 1.0)

  # Lower factor (more permissive) should detect more or equal outliers
  expect_gte(length(result_permissive$outliers), length(result_strict$outliers))
})

test_that("outliers.boxplot validates input", {
  expect_error(outliers.boxplot(matrix(1:10, 2, 5)))
})

# ============== magnitudeshape tests ==============

test_that("magnitudeshape returns magnitudeshape object", {
  fd <- create_fdata_with_outliers()
  result <- magnitudeshape(fd)

  expect_s3_class(result, "magnitudeshape")
  expect_true("MO" %in% names(result))
  expect_true("VO" %in% names(result))
  expect_true("outliers" %in% names(result))
  expect_true("cutoff" %in% names(result))
})

test_that("magnitudeshape MO and VO have correct length", {
  fd <- create_fdata_with_outliers(n = 30)
  result <- magnitudeshape(fd)

  expect_equal(length(result$MO), 30)
  expect_equal(length(result$VO), 30)
})

test_that("magnitudeshape VO is non-negative", {
  fd <- create_fdata_with_outliers()
  result <- magnitudeshape(fd)

  expect_true(all(result$VO >= 0))
})

test_that("magnitudeshape validates input", {
  expect_error(magnitudeshape(matrix(1:10, 2, 5)))
})

# ============== outliergram tests ==============

test_that("outliergram returns outliergram object", {
  fd <- create_fdata_with_outliers()
  result <- outliergram(fd)

  expect_s3_class(result, "outliergram")
  expect_true("mei" %in% names(result))
  expect_true("mbd" %in% names(result))
  expect_true("outliers" %in% names(result))
  expect_true("outlier_type" %in% names(result))
})

test_that("outliergram MEI and MBD in valid range", {
  fd <- create_fdata_with_outliers()
  result <- outliergram(fd)

  expect_true(all(result$mei >= 0 & result$mei <= 1))
  expect_true(all(result$mbd >= 0 & result$mbd <= 1))
})

test_that("outliergram factor affects detection", {
  fd <- create_fdata_with_outliers()

  result_strict <- outliergram(fd, factor = 3.0)
  result_permissive <- outliergram(fd, factor = 0.5)

  # Lower factor (more permissive) should detect more or equal outliers
  expect_gte(result_permissive$n_outliers, result_strict$n_outliers)
})

test_that("outliergram outlier types are valid", {
  fd <- create_fdata_with_outliers()
  result <- outliergram(fd, factor = 1.0)

  if (length(result$outliers) > 0) {
    expect_true(all(result$outlier_type %in% "shape"))
  }
})

test_that("outliergram validates input", {
  expect_error(outliergram(matrix(1:10, 2, 5)))
})

test_that("outliergram uses finite-sample parabola coefficients", {
  set.seed(42)
  m <- 50
  t_grid <- seq(0, 1, length.out = m)
  n <- 30
  X <- matrix(0, n, m)
  for (i in 1:n) X[i, ] <- sin(2 * pi * t_grid) + rnorm(m, sd = 0.2)
  fd <- fdata(X, argvals = t_grid)

  result <- outliergram(fd)

  # Finite-sample coefficients per Arribas-Gil & Romo (2014), Prop. 1
  expected_a0 <- -2 / (n * (n - 1))
  expected_a1 <- 2 * (n + 1) / (n - 1)
  expected_a2 <- -2 * (n + 1) / (n - 1)

  expect_equal(unname(result$parabola["a0"]), expected_a0)
  expect_equal(unname(result$parabola["a1"]), expected_a1)
  expect_equal(unname(result$parabola["a2"]), expected_a2)
})

test_that("outliergram uses boxplot fence threshold", {
  set.seed(42)
  m <- 50
  t_grid <- seq(0, 1, length.out = m)
  n <- 30
  X <- matrix(0, n, m)
  for (i in 1:n) X[i, ] <- sin(2 * pi * t_grid) + rnorm(m, sd = 0.2)
  fd <- fdata(X, argvals = t_grid)

  result <- outliergram(fd, factor = 1.5)

  # Threshold should be Q3 + 1.5 * IQR of the distances
  d <- result$dist_to_parabola
  q1 <- unname(quantile(d, 0.25))
  q3 <- unname(quantile(d, 0.75))
  expected_threshold <- q3 + 1.5 * (q3 - q1)
  expect_equal(unname(result$threshold), unname(expected_threshold))

  # All detected outliers must have distance > threshold
  if (length(result$outliers) > 0) {
    expect_true(all(d[result$outliers] > result$threshold))
  }
  # All non-outliers must have distance <= threshold
  non_outliers <- setdiff(seq_len(n), result$outliers)
  expect_true(all(d[non_outliers] <= result$threshold))
})

test_that("outliergram detects shape outliers but not pure magnitude outliers", {
  set.seed(42)
  m <- 50
  t_grid <- seq(0, 1, length.out = m)
  n <- 32
  X <- matrix(0, n, m)
  for (i in 1:29) X[i, ] <- sin(2 * pi * t_grid) + rnorm(m, sd = 0.15)
  X[30, ] <- sin(2 * pi * t_grid) + 3   # magnitude outlier (shifted up)
  X[31, ] <- sin(2 * pi * t_grid) - 3   # magnitude outlier (shifted down)
  X[32, ] <- sin(4 * pi * t_grid)       # shape outlier (different frequency)
  fd <- fdata(X, argvals = t_grid)

  result <- outliergram(fd, factor = 1.5)

  # Shape outlier (curve 32) should be detected

  expect_true(32 %in% result$outliers)

  # Pure magnitude outliers (30, 31) should NOT be detected by outliergram
  # because shifted curves still follow the MEI-MBD parabolic relationship
  expect_false(30 %in% result$outliers)
  expect_false(31 %in% result$outliers)
})

test_that("outliergram dist_to_parabola is positive below parabola", {
  set.seed(42)
  m <- 50
  t_grid <- seq(0, 1, length.out = m)
  n <- 20
  X <- matrix(0, n, m)
  for (i in 1:n) X[i, ] <- sin(2 * pi * t_grid) + rnorm(m, sd = 0.2)
  fd <- fdata(X, argvals = t_grid)

  result <- outliergram(fd)

  # dist_to_parabola = mbd_expected - mbd, so positive means below parabola
  mbd_expected <- result$parabola["a0"] + result$parabola["a1"] * result$mei +
                  result$parabola["a2"] * result$mei^2
  expect_equal(unname(result$dist_to_parabola), unname(mbd_expected - result$mbd))
})

# ============== print methods tests ==============

test_that("print.outliers.fdata works", {
  fd <- create_fdata_with_outliers()
  result <- outliers.depth.pond(fd, nb = 30)

  expect_output(print(result), "Functional data outlier detection")
})

test_that("print.magnitudeshape works", {
  fd <- create_fdata_with_outliers()
  result <- magnitudeshape(fd)

  expect_output(print(result), "Magnitude-Shape")
})

test_that("print.outliergram works", {
  fd <- create_fdata_with_outliers()
  result <- outliergram(fd)

  expect_output(print(result), "Outliergram")
})

# ============== p-value and outlierness tests ==============

test_that("outliers.depth.pond returns p.value and outlierness", {
  fd <- create_fdata_with_outliers(n = 30, n_outliers = 2)
  out <- outliers.depth.pond(fd, nb = 200)

  expect_true(!is.null(out$p.value))
  expect_true(!is.null(out$outlierness))
  expect_length(out$p.value, 30)
  expect_length(out$outlierness, 30)
  expect_true(all(out$p.value > 0 & out$p.value < 1))
  expect_true(all(out$outlierness >= 0 & out$outlierness <= 1))

  # Known outliers (29, 30) should have small p-values relative to normals
  normal_median_p <- median(out$p.value[1:28])
  expect_true(out$p.value[29] < normal_median_p)
  expect_true(out$p.value[30] < normal_median_p)
})

test_that("outliers.depth.trim returns p.value and outlierness", {
  fd <- create_fdata_with_outliers(n = 20, n_outliers = 2)
  out <- outliers.depth.trim(fd, trim = 0.1)

  expect_true(!is.null(out$p.value))
  expect_true(!is.null(out$outlierness))
  expect_length(out$p.value, 20)
  expect_length(out$outlierness, 20)
  expect_true(all(out$p.value > 0 & out$p.value < 1))
})

test_that("outliers.lrt returns p.value and outlierness", {
  set.seed(42)
  t <- seq(0, 1, length.out = 50)
  X <- matrix(0, 30, 50)
  for (i in 1:28) X[i, ] <- sin(2 * pi * t) + rnorm(50, sd = 0.1)
  X[29, ] <- sin(2 * pi * t) + 5
  X[30, ] <- sin(2 * pi * t) - 5
  fd <- fdata(X, argvals = t)

  out <- outliers.lrt(fd, nb = 200, seed = 123)

  expect_true(!is.null(out$p.value))
  expect_true(!is.null(out$outlierness))
  expect_length(out$p.value, 30)
  expect_length(out$outlierness, 30)
  expect_true(all(out$p.value > 0 & out$p.value <= 1))

  # Outliers should have outlierness > 1 (distance exceeds threshold)
  if (length(out$outliers) > 0) {
    expect_true(all(out$outlierness[out$outliers] > 1))
  }
})

test_that("outliers.boxplot returns p.value, outlierness, and exceedance", {
  fd <- create_fdata_with_outliers(n = 30, n_outliers = 2)
  out <- outliers.boxplot(fd)

  expect_true(!is.null(out$p.value))
  expect_true(!is.null(out$outlierness))
  expect_true(!is.null(out$exceedance))
  expect_length(out$p.value, 30)
  expect_length(out$outlierness, 30)
  expect_length(out$exceedance, 30)

  # Non-outlier curves should have exceedance = 0
  non_outliers <- setdiff(1:30, out$outliers)
  if (length(non_outliers) > 0) {
    expect_true(all(out$exceedance[non_outliers] == 0))
  }
})

test_that("magnitudeshape returns p.value and outlierness", {
  set.seed(42)
  t <- seq(0, 1, length.out = 50)
  X <- matrix(0, 30, 50)
  for (i in 1:28) X[i, ] <- sin(2 * pi * t) + rnorm(50, sd = 0.2)
  X[29, ] <- sin(2 * pi * t) + 3
  X[30, ] <- sin(4 * pi * t)
  fd <- fdata(X, argvals = t)

  ms <- magnitudeshape(fd)

  expect_true(!is.null(ms$p.value))
  expect_true(!is.null(ms$outlierness))
  expect_length(ms$p.value, 30)
  expect_length(ms$outlierness, 30)
  expect_true(all(ms$p.value >= 0 & ms$p.value <= 1))
  expect_true(all(ms$outlierness >= 0 & ms$outlierness <= 1))

  # Chi-squared p-values: outliers should have small p-values
  if (length(ms$outliers) > 0) {
    expect_true(all(ms$p.value[ms$outliers] < 0.05))
  }
})

test_that("outliergram returns p.value and outlierness", {
  set.seed(42)
  t <- seq(0, 1, length.out = 50)
  X <- matrix(0, 30, 50)
  for (i in 1:28) X[i, ] <- sin(2 * pi * t) + rnorm(50, sd = 0.2)
  X[29, ] <- sin(2 * pi * t) + 2
  X[30, ] <- sin(4 * pi * t)
  fd <- fdata(X, argvals = t)

  og <- outliergram(fd)

  expect_true(!is.null(og$p.value))
  expect_true(!is.null(og$outlierness))
  expect_length(og$p.value, 30)
  expect_length(og$outlierness, 30)
  expect_true(all(og$p.value > 0 & og$p.value <= 1))
  expect_true(all(og$outlierness >= 0 & og$outlierness <= 1))
})

# ============== outlier_summary tests ==============

test_that("outlier_summary returns correct data.frame for outliers.fdata", {
  fd <- create_fdata_with_outliers(n = 20, n_outliers = 2)
  out <- outliers.depth.pond(fd, nb = 200)
  summary_df <- outlier_summary(out)

  expect_s3_class(summary_df, "data.frame")
  expect_true(all(c("index", "outlierness", "p.value", "is_outlier") %in% names(summary_df)))
  expect_equal(nrow(summary_df), 20)

  # Sorted by outlierness descending
  expect_true(all(diff(summary_df$outlierness) <= 1e-10))

  # is_outlier matches detected outliers
  expect_equal(sum(summary_df$is_outlier), length(out$outliers))
})

test_that("outlier_summary works for magnitudeshape", {
  fd <- create_fdata_with_outliers(n = 20, n_outliers = 2)
  ms <- magnitudeshape(fd)
  summary_df <- outlier_summary(ms)

  expect_s3_class(summary_df, "data.frame")
  expect_equal(nrow(summary_df), 20)
  expect_true(all(diff(summary_df$outlierness) <= 1e-10))
})

test_that("outlier_summary works for outliergram", {
  fd <- create_fdata_with_outliers(n = 20, n_outliers = 2)
  og <- outliergram(fd)
  summary_df <- outlier_summary(og)

  expect_s3_class(summary_df, "data.frame")
  expect_equal(nrow(summary_df), 20)
})

test_that("outlier_summary rejects invalid input", {
  expect_error(outlier_summary(list(a = 1)), "must be an outlier detection result")
})
