# Tests for cv.fdata unified cross-validation framework

# =============================================================================
# Helper: create regression test data
# =============================================================================
create_regression_data <- function(n = 40, m = 30, seed = 42) {
  set.seed(seed)
  t_grid <- seq(0, 1, length.out = m)
  X <- matrix(0, n, m)
  for (i in seq_len(n)) {
    X[i, ] <- sin(2 * pi * t_grid * (1 + 0.2 * i / n)) + rnorm(m, sd = 0.1)
  }
  y <- rowMeans(X) + rnorm(n, sd = 0.1)
  list(fd = fdars::fdata(X, argvals = t_grid), y = y, n = n)
}

# =============================================================================
# Helper: create classification test data
# =============================================================================
create_classif_data_cv <- function(n_per_class = 15, m = 30, seed = 42) {
  set.seed(seed)
  t_grid <- seq(0, 1, length.out = m)
  data1 <- matrix(0, nrow = n_per_class, ncol = m)
  data2 <- matrix(0, nrow = n_per_class, ncol = m)
  for (i in seq_len(n_per_class)) {
    data1[i, ] <- sin(2 * pi * t_grid) + rnorm(m, 0, 0.15)
    data2[i, ] <- cos(2 * pi * t_grid) + rnorm(m, 0, 0.15)
  }
  data_mat <- rbind(data1, data2)
  y <- factor(rep(c("A", "B"), each = n_per_class))
  list(
    fd = fdars::fdata(data_mat, argvals = t_grid),
    y = y,
    n = 2 * n_per_class
  )
}

# =============================================================================
# Regression tests
# =============================================================================

test_that("cv.fdata regression produces correct structure", {
  d <- create_regression_data()
  result <- fdars::cv.fdata(d$fd, d$y,
    fit.fn = function(fd, y, ...) fdars::fregre.pc(fd, y, ncomp = 3),
    kfold = 5, seed = 123)

  expect_s3_class(result, "cv.fdata")
  expect_equal(result$type, "regression")
  expect_equal(result$kfold, 5)
  expect_length(result$oof.predictions, d$n)
  expect_length(result$folds, d$n)
  expect_length(result$fold.models, 5)
  expect_true(is.numeric(result$oof.predictions))

  # Metrics present
  expect_true("RMSE" %in% names(result$metrics))
  expect_true("MAE" %in% names(result$metrics))
  expect_true("R2" %in% names(result$metrics))
  expect_true(is.finite(result$metrics$RMSE))

  # Fold metrics
  expect_equal(nrow(result$fold.metrics), 5)
})

test_that("OOF predictions cover all observations exactly once", {
  d <- create_regression_data()
  result <- fdars::cv.fdata(d$fd, d$y,
    fit.fn = function(fd, y, ...) fdars::fregre.pc(fd, y, ncomp = 3),
    kfold = 5, seed = 123)

  # Every observation should have a prediction (no NAs)
  expect_true(all(!is.na(result$oof.predictions)))

  # Each observation assigned to exactly one fold
  expect_true(all(result$folds %in% 1:5))
  expect_equal(length(result$folds), d$n)

  # Each fold number appears at least once
  expect_equal(sort(unique(result$folds)), 1:5)
})

test_that("seed produces reproducible results", {
  d <- create_regression_data()
  fit_fn <- function(fd, y, ...) fdars::fregre.pc(fd, y, ncomp = 3)

  r1 <- fdars::cv.fdata(d$fd, d$y, fit.fn = fit_fn, kfold = 5, seed = 99)
  r2 <- fdars::cv.fdata(d$fd, d$y, fit.fn = fit_fn, kfold = 5, seed = 99)

  expect_identical(r1$folds, r2$folds)
  expect_equal(r1$oof.predictions, r2$oof.predictions)
  expect_equal(r1$metrics$RMSE, r2$metrics$RMSE)
})

# =============================================================================
# Classification tests
# =============================================================================

test_that("cv.fdata classification produces correct structure", {
  d <- create_classif_data_cv()
  result <- fdars::cv.fdata(d$fd, d$y,
    fit.fn = function(fd, y, ...) fdars::fregre.pc(fd, as.numeric(y), ncomp = 3),
    predict.fn = function(model, newdata) {
      preds <- predict(model, newdata)
      ifelse(preds < 1.5, "A", "B")
    },
    kfold = 5, type = "classification", seed = 42)

  expect_s3_class(result, "cv.fdata")
  expect_equal(result$type, "classification")
  expect_true(is.factor(result$oof.predictions))
  expect_equal(levels(result$oof.predictions), levels(d$y))
  expect_length(result$oof.predictions, d$n)
  expect_true("accuracy" %in% names(result$metrics))
  expect_true(!is.null(result$metrics$confusion))
})

test_that("auto-detect classification from factor y", {
  d <- create_classif_data_cv()
  result <- fdars::cv.fdata(d$fd, d$y,
    fit.fn = function(fd, y, ...) fdars::fregre.pc(fd, as.numeric(y), ncomp = 3),
    predict.fn = function(model, newdata) {
      preds <- predict(model, newdata)
      ifelse(preds < 1.5, "A", "B")
    },
    kfold = 3, seed = 42)

  # Should auto-detect classification since y is factor

  expect_equal(result$type, "classification")
})

# =============================================================================
# Stratification tests
# =============================================================================

test_that("stratified folds have balanced class distribution", {
  d <- create_classif_data_cv(n_per_class = 20)
  folds <- fdars::: .create_folds(d$y, kfold = 5, type = "classification",
                                   stratified = TRUE, seed = 42)

  # Each fold should have observations from both classes
  for (k in 1:5) {
    idx <- which(folds == k)
    classes_in_fold <- unique(d$y[idx])
    expect_true(length(classes_in_fold) == 2,
      info = paste("Fold", k, "missing a class"))
  }

  # Class proportions should be roughly equal across folds
  for (k in 1:5) {
    idx <- which(folds == k)
    prop_A <- mean(d$y[idx] == "A")
    expect_true(prop_A > 0.2 && prop_A < 0.8,
      info = paste("Fold", k, "has imbalanced classes"))
  }
})

test_that("stratified regression folds have balanced y distribution", {
  d <- create_regression_data(n = 50)
  folds <- fdars::: .create_folds(d$y, kfold = 5, type = "regression",
                                   stratified = TRUE, seed = 42)

  # Fold means of y should be roughly similar
  fold_means <- tapply(d$y, folds, mean)
  overall_mean <- mean(d$y)
  # Each fold mean should be within 1 SD of overall mean
  overall_sd <- sd(d$y)
  for (k in 1:5) {
    expect_true(abs(fold_means[k] - overall_mean) < 2 * overall_sd,
      info = paste("Fold", k, "has very different mean"))
  }
})

# =============================================================================
# Input validation tests
# =============================================================================

test_that("cv.fdata validates inputs", {
  d <- create_regression_data()

  expect_error(fdars::cv.fdata(d$fd$data, d$y, fit.fn = identity),
               "fdata")
  expect_error(fdars::cv.fdata(d$fd, d$y[1:5], fit.fn = identity),
               "must equal")
  expect_error(fdars::cv.fdata(d$fd, d$y, fit.fn = "not_a_function"),
               "must be a function")
  expect_error(fdars::cv.fdata(d$fd, d$y, fit.fn = identity, kfold = 1),
               "must be between")
})

# =============================================================================
# Print method test
# =============================================================================

test_that("print.cv.fdata works", {
  d <- create_regression_data()
  result <- fdars::cv.fdata(d$fd, d$y,
    fit.fn = function(fd, y, ...) fdars::fregre.pc(fd, y, ncomp = 3),
    kfold = 5, seed = 42)

  out <- capture.output(print(result))
  expect_true(any(grepl("K-Fold", out)))
  expect_true(any(grepl("RMSE", out)))
  expect_true(any(grepl("regression", out)))
})

# =============================================================================
# Repeated Cross-Validation tests (nrep > 1)
# =============================================================================

test_that("nrep = 1 produces identical results to default (backward compat)", {
  d <- create_regression_data()
  fit_fn <- function(fd, y, ...) fdars::fregre.pc(fd, y, ncomp = 3)

  r_default <- fdars::cv.fdata(d$fd, d$y, fit.fn = fit_fn, kfold = 5,
                                seed = 42)
  r_nrep1 <- fdars::cv.fdata(d$fd, d$y, fit.fn = fit_fn, kfold = 5,
                               nrep = 1, seed = 42)

  expect_equal(r_default$oof.predictions, r_nrep1$oof.predictions)
  expect_identical(r_default$folds, r_nrep1$folds)
  expect_equal(r_default$metrics$RMSE, r_nrep1$metrics$RMSE)
  expect_null(r_nrep1$oof.matrix)
  expect_null(r_nrep1$oof.sd)
  expect_null(r_nrep1$nrep)
})

test_that("nrep = 3 returns oof.matrix with correct dimensions", {
  d <- create_regression_data()
  result <- fdars::cv.fdata(d$fd, d$y,
    fit.fn = function(fd, y, ...) fdars::fregre.pc(fd, y, ncomp = 3),
    kfold = 5, nrep = 3, seed = 42)

  expect_equal(nrow(result$oof.matrix), d$n)
  expect_equal(ncol(result$oof.matrix), 3)
  expect_true(is.matrix(result$oof.matrix))
})

test_that("oof.sd is numeric(n) with no NAs for repeated CV", {
  d <- create_regression_data()
  result <- fdars::cv.fdata(d$fd, d$y,
    fit.fn = function(fd, y, ...) fdars::fregre.pc(fd, y, ncomp = 3),
    kfold = 5, nrep = 3, seed = 42)

  expect_true(is.numeric(result$oof.sd))
  expect_length(result$oof.sd, d$n)
  expect_true(all(!is.na(result$oof.sd)))
  expect_true(all(result$oof.sd >= 0))
})

test_that("folds.matrix columns differ across reps", {
  d <- create_regression_data()
  result <- fdars::cv.fdata(d$fd, d$y,
    fit.fn = function(fd, y, ...) fdars::fregre.pc(fd, y, ncomp = 3),
    kfold = 5, nrep = 3, seed = 42)

  expect_equal(dim(result$folds.matrix), c(d$n, 3))
  # Fold assignments should differ between reps
  expect_false(identical(result$folds.matrix[, 1], result$folds.matrix[, 2]))
  expect_false(identical(result$folds.matrix[, 1], result$folds.matrix[, 3]))
})

test_that("rep.metrics has nrep rows", {
  d <- create_regression_data()
  result <- fdars::cv.fdata(d$fd, d$y,
    fit.fn = function(fd, y, ...) fdars::fregre.pc(fd, y, ncomp = 3),
    kfold = 5, nrep = 3, seed = 42)

  expect_equal(nrow(result$rep.metrics), 3)
  expect_true("RMSE" %in% names(result$rep.metrics))
  expect_true("MAE" %in% names(result$rep.metrics))
  expect_true("R2" %in% names(result$rep.metrics))
  expect_true(all(result$rep.metrics$rep == 1:3))
})

test_that("metrics.summary has mean and sd", {
  d <- create_regression_data()
  result <- fdars::cv.fdata(d$fd, d$y,
    fit.fn = function(fd, y, ...) fdars::fregre.pc(fd, y, ncomp = 3),
    kfold = 5, nrep = 3, seed = 42)

  ms <- result$metrics.summary
  expect_true(is.list(ms))
  expect_true(all(c("RMSE", "MAE", "R2") %in% names(ms)))
  for (metric in c("RMSE", "MAE", "R2")) {
    expect_true("mean" %in% names(ms[[metric]]))
    expect_true("sd" %in% names(ms[[metric]]))
    expect_true(is.finite(ms[[metric]]["mean"]))
    expect_true(is.finite(ms[[metric]]["sd"]))
  }
})

test_that("seed reproducibility for repeated CV", {
  d <- create_regression_data()
  fit_fn <- function(fd, y, ...) fdars::fregre.pc(fd, y, ncomp = 3)

  r1 <- fdars::cv.fdata(d$fd, d$y, fit.fn = fit_fn, kfold = 5, nrep = 3,
                          seed = 99)
  r2 <- fdars::cv.fdata(d$fd, d$y, fit.fn = fit_fn, kfold = 5, nrep = 3,
                          seed = 99)

  expect_equal(r1$oof.matrix, r2$oof.matrix)
  expect_identical(r1$folds.matrix, r2$folds.matrix)
  expect_equal(r1$oof.sd, r2$oof.sd)
  expect_equal(r1$metrics$RMSE, r2$metrics$RMSE)
})

test_that("fold.models is list of nrep lists for repeated CV", {
  d <- create_regression_data()
  result <- fdars::cv.fdata(d$fd, d$y,
    fit.fn = function(fd, y, ...) fdars::fregre.pc(fd, y, ncomp = 3),
    kfold = 5, nrep = 3, seed = 42)

  expect_true(is.list(result$fold.models))
  expect_length(result$fold.models, 3)
  for (r in 1:3) {
    expect_true(is.list(result$fold.models[[r]]))
    expect_length(result$fold.models[[r]], 5)
  }
})

test_that("print includes 'repetitions' when nrep > 1", {
  d <- create_regression_data()
  result <- fdars::cv.fdata(d$fd, d$y,
    fit.fn = function(fd, y, ...) fdars::fregre.pc(fd, y, ncomp = 3),
    kfold = 5, nrep = 3, seed = 42)

  out <- capture.output(print(result))
  expect_true(any(grepl("Repeated", out)))
  expect_true(any(grepl("repetitions", out)))
  expect_true(any(grepl("Per-repetition", out)))
  expect_true(any(grepl("variability", out, ignore.case = TRUE)))
})

test_that("classification with nrep > 1 uses majority vote", {
  d <- create_classif_data_cv()
  result <- fdars::cv.fdata(d$fd, d$y,
    fit.fn = function(fd, y, ...) fdars::fregre.pc(fd, as.numeric(y),
                                                     ncomp = 3),
    predict.fn = function(model, newdata) {
      preds <- predict(model, newdata)
      ifelse(preds < 1.5, "A", "B")
    },
    kfold = 5, nrep = 3, type = "classification", seed = 42)

  expect_equal(result$nrep, 3)
  expect_true(is.factor(result$oof.predictions))
  expect_equal(levels(result$oof.predictions), levels(d$y))
  expect_length(result$oof.predictions, d$n)

  # oof.sd should be disagreement proportion (0 to 1)
  expect_true(all(result$oof.sd >= 0 & result$oof.sd <= 1))

  # oof.matrix should be character matrix
  expect_equal(dim(result$oof.matrix), c(d$n, 3))

  # rep.metrics should have accuracy
  expect_equal(nrow(result$rep.metrics), 3)
  expect_true("accuracy" %in% names(result$rep.metrics))
})
