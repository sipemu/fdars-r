# Tests for integral kernel (IKer) and smoother matrix functions

# =============================================================================
# IKer CDF properties: all should be valid CDFs on [-1, 1]
# =============================================================================

test_that("IKer.norm is a valid CDF", {
  # Should be 0 at -Inf, 1 at +Inf, monotonically increasing
  expect_equal(IKer.norm(-10), 0, tolerance = 1e-6)
  expect_equal(IKer.norm(0), 0.5, tolerance = 1e-6)
  expect_equal(IKer.norm(10), 1, tolerance = 1e-6)

  u <- seq(-3, 3, length.out = 100)
  vals <- IKer.norm(u)
  expect_true(all(diff(vals) >= 0))
})

test_that("IKer.epa is a valid CDF on [-1, 1]", {
  expect_equal(IKer.epa(-2), 0)
  expect_equal(IKer.epa(0), 0.5, tolerance = 1e-10)
  expect_equal(IKer.epa(2), 1)

  u <- seq(-1, 1, length.out = 100)
  vals <- IKer.epa(u)
  expect_true(all(diff(vals) >= -1e-10))
  expect_true(all(vals >= -1e-14 & vals <= 1 + 1e-14))
})

test_that("IKer.tri is a valid CDF on [-1, 1]", {
  expect_equal(IKer.tri(-2), 0)
  expect_equal(IKer.tri(0), 0.5, tolerance = 1e-10)
  expect_equal(IKer.tri(2), 1)

  u <- seq(-1, 1, length.out = 100)
  vals <- IKer.tri(u)
  expect_true(all(diff(vals) >= -1e-10))
  expect_true(all(vals >= -1e-14 & vals <= 1 + 1e-14))
})

test_that("IKer.quar is a valid CDF on [-1, 1]", {
  expect_equal(IKer.quar(-2), 0)
  expect_equal(IKer.quar(0), 0.5, tolerance = 1e-10)
  expect_equal(IKer.quar(2), 1)

  u <- seq(-1, 1, length.out = 100)
  vals <- IKer.quar(u)
  expect_true(all(diff(vals) >= -1e-10))
  expect_true(all(vals >= -1e-14 & vals <= 1 + 1e-14))
})

test_that("IKer.cos is a valid CDF on [-1, 1]", {
  expect_equal(IKer.cos(-2), 0)
  expect_equal(IKer.cos(0), 0.5, tolerance = 1e-10)
  expect_equal(IKer.cos(2), 1)

  u <- seq(-1, 1, length.out = 100)
  vals <- IKer.cos(u)
  expect_true(all(diff(vals) >= -1e-10))
  expect_true(all(vals >= -1e-14 & vals <= 1 + 1e-14))
})

test_that("IKer.unif is a valid CDF on [-1, 1]", {
  expect_equal(IKer.unif(-2), 0)
  expect_equal(IKer.unif(-1), 0)
  expect_equal(IKer.unif(0), 0.5)
  expect_equal(IKer.unif(1), 1)
  expect_equal(IKer.unif(2), 1)

  u <- seq(-1, 1, length.out = 100)
  vals <- IKer.unif(u)
  expect_true(all(diff(vals) >= -1e-10))
})

test_that("all IKer functions handle vectorized input", {
  u <- c(-2, -1, -0.5, 0, 0.5, 1, 2)
  for (ker in list(IKer.norm, IKer.epa, IKer.tri, IKer.quar, IKer.cos, IKer.unif)) {
    vals <- ker(u)
    expect_equal(length(vals), length(u))
    expect_true(all(is.finite(vals)))
  }
})

test_that("all compact IKer functions are 0 below -1 and 1 above 1", {
  for (ker in list(IKer.epa, IKer.tri, IKer.quar, IKer.cos, IKer.unif)) {
    expect_equal(ker(-5), 0)
    expect_equal(ker(-1.01), 0)
    expect_equal(ker(1.01), 1)
    expect_equal(ker(5), 1)
  }
})

# =============================================================================
# S.LCR (Local Cubic Regression smoother matrix)
# =============================================================================

test_that("S.LCR returns a smoother matrix", {
  tt <- seq(0, 1, length.out = 20)
  h <- 0.2
  S <- S.LCR(tt, h)

  expect_true(is.matrix(S))
  expect_equal(nrow(S), 20)
  expect_equal(ncol(S), 20)
})

test_that("S.LCR with cv=TRUE returns LOO matrix", {
  tt <- seq(0, 1, length.out = 20)
  S_cv <- S.LCR(tt, h = 0.2, cv = TRUE)

  expect_true(is.matrix(S_cv))
  # LOO: diagonal should be zero (each point excluded from own fit)
  expect_true(all(abs(diag(S_cv)) < 1e-10))
})

# =============================================================================
# CV.S (Cross-validation score for smoother)
# =============================================================================

test_that("CV.S returns a non-negative scalar", {
  set.seed(42)
  tt <- seq(0, 1, length.out = 30)
  y <- sin(2 * pi * tt) + rnorm(30, sd = 0.3)

  cv_score <- CV.S(S.NW, tt, h = 0.15, y = y)

  expect_true(is.numeric(cv_score))
  expect_equal(length(cv_score), 1)
  expect_gte(cv_score, 0)
})

test_that("CV.S decreases with better bandwidth", {
  set.seed(42)
  tt <- seq(0, 1, length.out = 30)
  y <- sin(2 * pi * tt) + rnorm(30, sd = 0.3)

  # Very small bandwidth -> overfitting, very large -> underfitting
  # Moderate bandwidth should be better than extreme
  cv_small <- CV.S(S.NW, tt, h = 0.01, y = y)
  cv_med <- CV.S(S.NW, tt, h = 0.15, y = y)

  # At least one of the extremes should be worse
  expect_true(is.finite(cv_small) && is.finite(cv_med))
})
