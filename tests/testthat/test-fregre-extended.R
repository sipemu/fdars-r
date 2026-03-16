# Extended tests for fregre.R: scalar covariates, CV, prediction, rp.stat

# =============================================================================
# Helper
# =============================================================================

.test_fregre_data <- function(n = 60, m = 20, seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  z <- matrix(rnorm(n * 2), n, 2)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(m, sd = 0.1)
    y[i] <- 2 * amp + 0.5 * z[i, 1] + rnorm(1, sd = 0.3)
  }
  fd <- fdata(X, argvals = t)
  list(fd = fd, y = y, z = z)
}

# =============================================================================
# fregre.lm with scalar covariates
# =============================================================================

test_that("fregre.lm with scalar covariates returns correct structure", {
  d <- .test_fregre_data()
  fit <- fregre.lm(d$fd, d$y, scalar.covariates = d$z, ncomp = 3)

  expect_s3_class(fit, "fregre.lm")
  expect_equal(length(fit$fitted.values), 60)
  expect_equal(length(fit$residuals), 60)
  expect_true(is.numeric(fit$r.squared))
  expect_true(fit$r.squared >= 0 && fit$r.squared <= 1)
  expect_true(length(fit$gamma) > 0)  # scalar covariate coefficients
})

test_that("fregre.lm fitted + residuals = y", {
  d <- .test_fregre_data()
  fit <- fregre.lm(d$fd, d$y, ncomp = 3)

  expect_equal(fit$fitted.values + fit$residuals, d$y, tolerance = 1e-10)
})

test_that("fregre.lm stores model diagnostics", {
  d <- .test_fregre_data()
  fit <- fregre.lm(d$fd, d$y, ncomp = 3)

  expect_true(is.finite(fit$aic))
  expect_true(is.finite(fit$bic))
  expect_true(is.finite(fit$gcv))
  expect_true(fit$residual.se > 0)
  expect_true(is.numeric(fit$r.squared.adj))
  expect_true(length(fit$std.errors) > 0)
})

test_that("fregre.lm auto-selects ncomp when NULL", {
  d <- .test_fregre_data()
  fit <- fregre.lm(d$fd, d$y)

  expect_s3_class(fit, "fregre.lm")
  expect_true(fit$ncomp > 0)
})

test_that("fregre.lm validates input", {
  expect_error(fregre.lm("not_fdata", 1:10), "fdata")

  d <- .test_fregre_data()
  expect_error(fregre.lm(d$fd, 1:5), "Length of y")
})

# =============================================================================
# predict.fregre.lm
# =============================================================================

test_that("predict.fregre.lm returns fitted when newdata=NULL", {
  d <- .test_fregre_data()
  fit <- fregre.lm(d$fd, d$y, ncomp = 3)

  pred <- predict(fit)
  expect_equal(pred, fit$fitted.values)
})

test_that("predict.fregre.lm with new functional data", {
  d <- .test_fregre_data()
  fit <- fregre.lm(d$fd, d$y, ncomp = 3)

  pred <- predict(fit, d$fd[1:10, ])
  expect_equal(length(pred), 10)
  expect_true(all(is.finite(pred)))
})

test_that("predict.fregre.lm with scalar covariates", {
  d <- .test_fregre_data()
  fit <- fregre.lm(d$fd, d$y, scalar.covariates = d$z, ncomp = 3)

  pred <- predict(fit, d$fd[1:10, ], new.scalar = d$z[1:10, ])
  expect_equal(length(pred), 10)
  expect_true(all(is.finite(pred)))
})

test_that("predict.fregre.lm predictions are reasonable", {
  d <- .test_fregre_data()
  fit <- fregre.lm(d$fd, d$y, ncomp = 3)

  pred <- predict(fit, d$fd)
  # Predictions should correlate well with y for signal data
  r <- cor(pred, d$y)
  expect_gt(r, 0.5)
})

# =============================================================================
# fregre.lm.cv
# =============================================================================

test_that("fregre.lm.cv selects optimal ncomp", {
  d <- .test_fregre_data()
  cv <- fregre.lm.cv(d$fd, d$y, k.range = 1:5, nfold = 5)

  expect_true(!is.null(cv$optimal.k))
  expect_true(cv$optimal.k >= 1 && cv$optimal.k <= 5)
  expect_true(!is.null(cv$cv.errors))
  expect_s3_class(cv$model, "fregre.lm")
})

test_that("fregre.lm.cv validates input", {
  expect_error(fregre.lm.cv("not_fdata", 1:10), "fdata")
})

# =============================================================================
# fregre.pc
# =============================================================================

test_that("fregre.pc returns correct structure", {
  d <- .test_fregre_data()
  fit <- fregre.pc(d$fd, d$y, ncomp = 3)

  expect_s3_class(fit, "fregre.fd")
  expect_equal(length(fit$fitted.values), 60)
  expect_equal(length(fit$residuals), 60)
  expect_equal(fit$ncomp, 3)
  expect_true(is.numeric(fit$intercept))
})

test_that("fregre.pc predict works", {
  d <- .test_fregre_data()
  fit <- fregre.pc(d$fd, d$y, ncomp = 3)

  pred <- predict(fit, d$fd[1:10, ])
  expect_equal(length(pred), 10)
  expect_true(all(is.finite(pred)))
})

# =============================================================================
# fregre.basis
# =============================================================================

test_that("fregre.basis returns correct structure", {
  d <- .test_fregre_data()
  fit <- fregre.basis(d$fd, d$y)

  expect_s3_class(fit, "fregre.fd")
  expect_equal(length(fit$fitted.values), 60)
  expect_true(is.numeric(fit$r.squared))
})

test_that("fregre.basis predict works", {
  d <- .test_fregre_data()
  fit <- fregre.basis(d$fd, d$y)

  pred <- predict(fit, d$fd[1:10, ])
  expect_equal(length(pred), 10)
  expect_true(all(is.finite(pred)))
})

# =============================================================================
# fregre.np
# =============================================================================

test_that("fregre.np returns correct structure", {
  d <- .test_fregre_data()
  fit <- fregre.np(d$fd, d$y)

  expect_s3_class(fit, "fregre.np")
  expect_equal(length(fit$fitted.values), 60)
  expect_true(is.numeric(fit$fitted.values))
  expect_true(fit$h.opt > 0)
})

test_that("fregre.np predict with new data", {
  d <- .test_fregre_data()
  fit <- fregre.np(d$fd, d$y)

  pred <- predict(fit, d$fd[1:10, ])
  expect_equal(length(pred), 10)
  expect_true(all(is.finite(pred)))
})

test_that("fregre.np with knn", {
  d <- .test_fregre_data()
  fit <- fregre.np(d$fd, d$y, knn = 10)

  expect_s3_class(fit, "fregre.np")
  expect_equal(length(fit$fitted.values), 60)
})

# =============================================================================
# fregre.pc.cv / fregre.basis.cv / fregre.np.cv
# =============================================================================

test_that("fregre.pc.cv returns correct structure", {
  d <- .test_fregre_data()
  cv <- fregre.pc.cv(d$fd, d$y, kfold = 3, ncomp.range = 1:4, seed = 42)

  expect_true(!is.null(cv$optimal.ncomp))
  expect_true(cv$optimal.ncomp >= 1 && cv$optimal.ncomp <= 4)
  expect_true(!is.null(cv$cv.errors))
  expect_s3_class(cv$model, "fregre.fd")
})

test_that("fregre.basis.cv returns correct structure", {
  d <- .test_fregre_data()
  lambdas <- c(0.01, 0.1, 1)
  cv <- fregre.basis.cv(d$fd, d$y, kfold = 3, lambda.range = lambdas, seed = 42)

  expect_true(!is.null(cv$optimal.lambda))
  expect_true(!is.null(cv$cv.errors))
})

test_that("fregre.np.cv returns correct structure", {
  d <- .test_fregre_data()
  h_range <- c(0.05, 0.1, 0.5)
  cv <- fregre.np.cv(d$fd, d$y, kfold = 3, h.range = h_range, seed = 42)

  expect_true(!is.null(cv$optimal.h))
  expect_true(!is.null(cv$cv.errors))
})

# =============================================================================
# rp.stat
# =============================================================================

test_that("rp.stat returns CvM and KS statistics", {
  set.seed(42)
  n <- 30
  n_proj <- 5
  proj <- unlist(replicate(n_proj, sample(n), simplify = FALSE))
  resid <- rnorm(n)

  result <- rp.stat(proj, resid, n.proj = n_proj)

  expect_true(is.numeric(result$cvm))
  expect_true(is.numeric(result$ks))
  expect_equal(length(result$cvm), n_proj)
  expect_equal(length(result$ks), n_proj)
  expect_true(all(result$cvm >= 0))
  expect_true(all(result$ks >= 0))
})

# =============================================================================
# print.fregre.fd / print.fregre.np
# =============================================================================

test_that("print.fregre.fd works for fregre.pc", {
  d <- .test_fregre_data()
  fit <- fregre.pc(d$fd, d$y, ncomp = 3)
  expect_output(print(fit), "Functional")
})

test_that("print.fregre.np works", {
  d <- .test_fregre_data()
  fit <- fregre.np(d$fd, d$y)
  expect_output(print(fit), "Nonparametric")
})
