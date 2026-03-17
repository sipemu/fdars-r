# Tests for alignment Phases 6-7: lambda CV, warp statistics, phase boxplot,
# warp inverse, penalized alignment, alignment diagnostics

# Helper: create curves with phase variation
.test_phase_data <- function(n = 10, m = 30, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  for (i in 1:n) {
    phase <- runif(1, -0.1, 0.1)
    X[i, ] <- sin(2 * pi * (t + phase)) + rnorm(m, sd = 0.05)
  }
  fdata(X, argvals = t)
}

# =============================================================================
# elastic.lambda.cv
# =============================================================================

test_that("elastic.lambda.cv returns best_lambda and cv_scores", {
  fd <- .test_phase_data(n = 10, m = 30)
  lambdas <- 10^seq(-2, 1, length.out = 5)

  cv <- elastic.lambda.cv(fd, lambdas = lambdas, n.folds = 3,
                            max.iter = 5, seed = 42)

  expect_s3_class(cv, "lambda.cv")
  expect_true(is.numeric(cv$best.lambda))
  expect_true(cv$best.lambda > 0)
  expect_true(is.numeric(cv$cv.scores))
  expect_equal(length(cv$cv.scores), length(lambdas))
  expect_true(all(is.finite(cv$cv.scores)))
  expect_equal(length(cv$lambdas), length(lambdas))
})

test_that("elastic.lambda.cv best_lambda is one of the candidates", {
  fd <- .test_phase_data(n = 10, m = 30)
  lambdas <- c(0.01, 0.1, 1.0, 10.0)

  cv <- elastic.lambda.cv(fd, lambdas = lambdas, n.folds = 3,
                            max.iter = 5, seed = 42)

  expect_true(cv$best.lambda %in% cv$lambdas)
})

test_that("elastic.lambda.cv validates input", {
  expect_error(elastic.lambda.cv("not_fdata"), "fdata")
})

test_that("elastic.lambda.cv print method works", {
  fd <- .test_phase_data(n = 8, m = 25)
  cv <- elastic.lambda.cv(fd, lambdas = c(0.1, 1.0, 10.0),
                            n.folds = 2, max.iter = 3, seed = 1)
  expect_output(print(cv), "Lambda Cross-Validation")
})

# =============================================================================
# warp.statistics
# =============================================================================

test_that("warp.statistics returns all expected fields", {
  fd <- .test_phase_data(n = 10, m = 30)
  km <- karcher.mean(fd, max.iter = 5, tol = 1e-3)
  ws <- warp.statistics(km)

  expect_s3_class(ws, "warp.statistics")
  expect_true(is.numeric(ws$mean))
  expect_equal(length(ws$mean), 30)
  expect_true(is.numeric(ws$variance))
  expect_equal(length(ws$variance), 30)
  expect_true(is.numeric(ws$std.dev))
  expect_equal(length(ws$std.dev), 30)
  expect_true(is.numeric(ws$lower.band))
  expect_equal(length(ws$lower.band), 30)
  expect_true(is.numeric(ws$upper.band))
  expect_equal(length(ws$upper.band), 30)
  expect_true(is.numeric(ws$karcher.mean.warp))
  expect_true(is.numeric(ws$geodesic.distances))
  expect_equal(length(ws$geodesic.distances), 10)
  expect_true(all(ws$geodesic.distances >= 0))
})

test_that("warp.statistics bands are ordered correctly", {
  fd <- .test_phase_data(n = 10, m = 30)
  km <- karcher.mean(fd, max.iter = 5, tol = 1e-3)
  ws <- warp.statistics(km)

  # Lower band <= mean <= upper band
  expect_true(all(ws$lower.band <= ws$mean + 1e-10))
  expect_true(all(ws$mean <= ws$upper.band + 1e-10))
})

test_that("warp.statistics validates input class", {
  expect_error(warp.statistics("not_karcher"), "karcher.mean")
})

test_that("warp.statistics print method works", {
  fd <- .test_phase_data(n = 8, m = 25)
  km <- karcher.mean(fd, max.iter = 3, tol = 1e-3)
  ws <- warp.statistics(km)
  expect_output(print(ws), "Warping Function Statistics")
})

# =============================================================================
# phase.boxplot
# =============================================================================

test_that("phase.boxplot returns median, outlier_indices, etc", {
  fd <- .test_phase_data(n = 12, m = 30)
  km <- karcher.mean(fd, max.iter = 5, tol = 1e-3)
  pb <- phase.boxplot(km)

  expect_s3_class(pb, "phase.boxplot")
  expect_true(is.numeric(pb$median))
  expect_equal(length(pb$median), 30)
  expect_true(is.numeric(pb$median.index))
  expect_gte(pb$median.index, 1)
  expect_lte(pb$median.index, 12)
  expect_true(is.numeric(pb$central.lower))
  expect_true(is.numeric(pb$central.upper))
  expect_true(is.numeric(pb$whisker.lower))
  expect_true(is.numeric(pb$whisker.upper))
  expect_true(is.numeric(pb$depths))
  expect_equal(length(pb$depths), 12)
  expect_true(all(pb$depths >= 0))
  expect_true(is.numeric(pb$outlier.indices) || is.integer(pb$outlier.indices))
  expect_equal(pb$factor, 1.5)
})

test_that("phase.boxplot central region is within whiskers", {
  fd <- .test_phase_data(n = 12, m = 30)
  km <- karcher.mean(fd, max.iter = 5, tol = 1e-3)
  pb <- phase.boxplot(km)

  # Central region should be within whiskers
  expect_true(all(pb$whisker.lower <= pb$central.lower + 1e-10))
  expect_true(all(pb$central.upper <= pb$whisker.upper + 1e-10))
})

test_that("phase.boxplot with different factor", {
  fd <- .test_phase_data(n = 10, m = 25)
  km <- karcher.mean(fd, max.iter = 3, tol = 1e-3)

  pb_small <- phase.boxplot(km, factor = 1.0)
  pb_large <- phase.boxplot(km, factor = 3.0)

  expect_equal(pb_small$factor, 1.0)
  expect_equal(pb_large$factor, 3.0)
})

test_that("phase.boxplot validates input class", {
  expect_error(phase.boxplot("not_karcher"), "karcher.mean")
})

test_that("phase.boxplot print method works", {
  fd <- .test_phase_data(n = 8, m = 25)
  km <- karcher.mean(fd, max.iter = 3, tol = 1e-3)
  pb <- phase.boxplot(km)
  expect_output(print(pb), "Phase Boxplot")
})

# =============================================================================
# warp.inverse
# =============================================================================

test_that("warp.inverse of identity warp is identity", {
  t <- seq(0, 1, length.out = 50)
  gamma <- t  # identity warp
  gamma_inv <- warp.inverse(gamma, t)

  expect_true(is.numeric(gamma_inv))
  expect_equal(length(gamma_inv), 50)
  expect_equal(gamma_inv, t, tolerance = 1e-6)
})

test_that("warp.inverse produces valid warp", {
  t <- seq(0, 1, length.out = 50)
  gamma <- t^2  # non-trivial warp

  gamma_inv <- warp.inverse(gamma, t)

  expect_true(is.numeric(gamma_inv))
  expect_equal(length(gamma_inv), 50)
  # Boundary conditions: gamma_inv(0) ~ 0, gamma_inv(1) ~ 1
  expect_equal(gamma_inv[1], 0, tolerance = 1e-4)
  expect_equal(gamma_inv[50], 1, tolerance = 1e-4)
  # Should be monotonically increasing
  diffs <- diff(gamma_inv)
  expect_true(all(diffs >= -1e-10))
})

test_that("warp.inverse is approximate left-inverse", {
  t <- seq(0, 1, length.out = 50)
  gamma <- t^2

  gamma_inv <- warp.inverse(gamma, t)

  # gamma(gamma_inv(t)) should be close to t
  # We approximate by interpolation
  composed <- approx(t, gamma, xout = gamma_inv)$y
  expect_equal(composed, t, tolerance = 0.05)
})

# =============================================================================
# warp.inverse.error
# =============================================================================

test_that("warp.inverse.error of identity warp is near zero", {
  t <- seq(0, 1, length.out = 50)
  gamma <- t  # identity warp

  err <- warp.inverse.error(gamma, t)

  expect_true(is.numeric(err))
  expect_equal(length(err), 1)
  expect_lt(err, 1e-6)
})

test_that("warp.inverse.error of non-trivial warp is small", {
  t <- seq(0, 1, length.out = 50)
  gamma <- t^2

  err <- warp.inverse.error(gamma, t)

  expect_true(is.numeric(err))
  expect_lt(err, 0.1)  # should be small for a smooth warp
})

test_that("warp.inverse.error is non-negative", {
  t <- seq(0, 1, length.out = 30)
  gamma <- t^1.5

  err <- warp.inverse.error(gamma, t)
  expect_gte(err, 0)
})

# =============================================================================
# elastic.pair.penalized
# =============================================================================

test_that("elastic.pair.penalized returns gamma, f.aligned, distance", {
  t <- seq(0, 1, length.out = 30)
  f1 <- sin(2 * pi * t)
  f2 <- sin(2 * pi * t^1.3)

  result <- elastic.pair.penalized(f1, f2, argvals = t)

  expect_true(is.list(result))
  expect_true(is.numeric(result$gamma))
  expect_equal(length(result$gamma), 30)
  expect_true(is.numeric(result$f.aligned))
  expect_equal(length(result$f.aligned), 30)
  expect_true(is.numeric(result$distance))
  expect_gte(result$distance, 0)
})

test_that("elastic.pair.penalized gamma is a valid warp", {
  t <- seq(0, 1, length.out = 30)
  f1 <- sin(2 * pi * t)
  f2 <- sin(2 * pi * t^1.3)

  result <- elastic.pair.penalized(f1, f2, argvals = t)

  # Boundary conditions
  expect_equal(result$gamma[1], 0, tolerance = 1e-4)
  expect_equal(result$gamma[30], 1, tolerance = 1e-4)
  # Monotonicity
  diffs <- diff(result$gamma)
  expect_true(all(diffs >= -1e-10))
})

test_that("elastic.pair.penalized identical curves have zero distance", {
  t <- seq(0, 1, length.out = 30)
  f <- sin(2 * pi * t)

  result <- elastic.pair.penalized(f, f, argvals = t)
  expect_lt(result$distance, 1e-4)
})

test_that("elastic.pair.penalized with different penalty types", {
  t <- seq(0, 1, length.out = 30)
  f1 <- sin(2 * pi * t)
  f2 <- sin(2 * pi * t^1.3)

  r1 <- elastic.pair.penalized(f1, f2, argvals = t, penalty = "first_order")
  r2 <- elastic.pair.penalized(f1, f2, argvals = t, penalty = "second_order")
  r3 <- elastic.pair.penalized(f1, f2, argvals = t, penalty = "combined")

  expect_true(is.list(r1))
  expect_true(is.list(r2))
  expect_true(is.list(r3))
  expect_true(r1$distance >= 0)
  expect_true(r2$distance >= 0)
  expect_true(r3$distance >= 0)
})

test_that("elastic.pair.penalized default argvals when NULL", {
  f1 <- sin(2 * pi * seq(0, 1, length.out = 30))
  f2 <- cos(2 * pi * seq(0, 1, length.out = 30))

  result <- elastic.pair.penalized(f1, f2)
  expect_true(is.list(result))
  expect_equal(length(result$gamma), 30)
})

# =============================================================================
# alignment.diagnostics
# =============================================================================

test_that("alignment.diagnostics returns health_score in [0,1]", {
  fd <- .test_phase_data(n = 10, m = 30)
  km <- karcher.mean(fd, max.iter = 5, tol = 1e-3)
  diag <- alignment.diagnostics(fd, km)

  expect_s3_class(diag, "alignment.diagnostics")
  expect_true(is.numeric(diag$health.score))
  expect_gte(diag$health.score, 0)
  expect_lte(diag$health.score, 1)
})

test_that("alignment.diagnostics returns per-curve diagnostics", {
  fd <- .test_phase_data(n = 10, m = 30)
  km <- karcher.mean(fd, max.iter = 5, tol = 1e-3)
  diag <- alignment.diagnostics(fd, km)

  expect_true(is.list(diag$diagnostics))
  expect_equal(length(diag$diagnostics), 10)

  # Check first diagnostic record
  d1 <- diag$diagnostics[[1]]
  expect_true(!is.null(d1$curve_index) || !is.null(d1$curve.index))
  expect_true(!is.null(d1$warp_complexity) || !is.null(d1$warp.complexity))
})

test_that("alignment.diagnostics flagged.indices is valid", {
  fd <- .test_phase_data(n = 10, m = 30)
  km <- karcher.mean(fd, max.iter = 5, tol = 1e-3)
  diag <- alignment.diagnostics(fd, km)

  expect_true(is.numeric(diag$flagged.indices) ||
              is.integer(diag$flagged.indices))
  expect_true(is.numeric(diag$n.flagged) || is.integer(diag$n.flagged))
  expect_gte(diag$n.flagged, 0)
  expect_lte(diag$n.flagged, 10)
  expect_equal(diag$n.flagged, length(diag$flagged.indices))
})

test_that("alignment.diagnostics validates input", {
  fd <- .test_phase_data(n = 10, m = 30)
  km <- karcher.mean(fd, max.iter = 3, tol = 1e-3)

  expect_error(alignment.diagnostics("not_fdata", km), "fdata")
  expect_error(alignment.diagnostics(fd, "not_karcher"), "karcher.mean")
})

test_that("alignment.diagnostics print method works", {
  fd <- .test_phase_data(n = 8, m = 25)
  km <- karcher.mean(fd, max.iter = 3, tol = 1e-3)
  diag <- alignment.diagnostics(fd, km)
  expect_output(print(diag), "Alignment Diagnostics")
})
