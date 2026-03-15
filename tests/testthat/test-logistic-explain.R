# Tests for logistic regression explainability functions

# Helper: create a fitted fregre.logistic model for testing
.test_logistic_model <- function(seed = 42) {
  set.seed(seed)
  t_grid <- seq(0, 1, length.out = 30)
  n <- 40
  X <- matrix(0, n, 30)
  for (i in 1:n) {
    amp <- if (i <= 20) rnorm(1, 1, 0.2) else rnorm(1, -1, 0.2)
    X[i, ] <- amp * sin(2 * pi * t_grid) + rnorm(30, sd = 0.1)
  }
  fd <- fdata(X, argvals = t_grid)
  y <- c(rep(0, 20), rep(1, 20))
  model <- functional.logistic(fd, y, ncomp = 3)
  list(model = model, data = fd, y = y)
}

# =============================================================================
# functional.logistic basics
# =============================================================================

test_that("functional.logistic returns correct structure", {
  env <- .test_logistic_model()
  m <- env$model

  expect_s3_class(m, "fregre.logistic")
  expect_true(is.numeric(m$intercept))
  expect_s3_class(m$beta.t, "fdata")
  expect_equal(length(m$probabilities), 40)
  expect_equal(length(m$predicted.classes), 40)
  expect_true(all(m$probabilities >= 0 & m$probabilities <= 1))
  expect_true(all(m$predicted.classes %in% c(0, 1)))
  expect_true(m$accuracy >= 0 && m$accuracy <= 1)
  expect_true(is.numeric(m$log.likelihood))
  expect_true(is.numeric(m$aic))
  expect_true(is.numeric(m$bic))
  expect_equal(m$ncomp, 3L)
})

test_that("predict.fregre.logistic works", {
  env <- .test_logistic_model()
  probs <- predict(env$model, env$data, type = "response")
  classes <- predict(env$model, env$data, type = "class")

  expect_equal(length(probs), 40)
  expect_true(all(probs >= 0 & probs <= 1))
  expect_equal(length(classes), 40)
  expect_true(all(classes %in% c(0, 1)))
})

test_that("predict.fregre.logistic with NULL newdata returns fitted", {
  env <- .test_logistic_model()
  probs <- predict(env$model, type = "response")
  expect_equal(probs, env$model$probabilities)
})

# =============================================================================
# Explainability methods with logistic model
# =============================================================================

test_that("fregre.beta.decomp works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.beta.decomp(env$model)
  expect_true(!is.null(result))
})

test_that("fregre.pointwise works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.pointwise(env$model)
  expect_true(!is.null(result))
  expect_true("importance" %in% names(result) || "importance_normalized" %in% names(result))
})

test_that("fregre.shap works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.shap(env$model, env$data)
  expect_true(!is.null(result))
  expect_true("values" %in% names(result) || "base_value" %in% names(result))
})

test_that("fregre.influence rejects logistic model", {
  env <- .test_logistic_model()
  expect_error(fregre.influence(env$model, env$data), "fregre\\.lm|fregre\\.robust")
})

test_that("fregre.vif works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.vif(env$model, env$data)
  expect_true(!is.null(result))
  expect_true("vif" %in% names(result))
})

test_that("fregre.dfbetas rejects logistic model", {
  env <- .test_logistic_model()
  expect_error(fregre.dfbetas(env$model, env$data), "fregre\\.lm|fregre\\.robust")
})

test_that("fregre.saliency works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.saliency(env$model, env$data)
  expect_true(!is.null(result))
  expect_true("saliency_map" %in% names(result) || "mean_absolute_saliency" %in% names(result))
})

test_that("fregre.pdp works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.pdp(env$model, env$data, component = 1, n.grid = 10)
  expect_true(!is.null(result))
  expect_true("grid_values" %in% names(result) || "pdp_curve" %in% names(result))
})

test_that("fregre.ale works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.ale(env$model, env$data, component = 1, n.bins = 10)
  expect_true(!is.null(result))
  expect_true("ale_values" %in% names(result) || "bin_midpoints" %in% names(result))
})

test_that("fregre.importance works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.importance(env$model, env$data, env$y, n.perm = 10, seed = 42)
  expect_true(!is.null(result))
  expect_true("importance" %in% names(result) || "baseline_metric" %in% names(result))
})

test_that("fregre.loo rejects logistic model", {
  env <- .test_logistic_model()
  expect_error(fregre.loo(env$model, env$data, env$y), "fregre\\.lm|fregre\\.robust")
})

test_that("fregre.sobol works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.sobol(env$model, env$data, env$y)
  expect_true(!is.null(result))
  expect_true("first_order" %in% names(result) || "total_order" %in% names(result))
})

test_that("fregre.domain works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.domain(env$model, window.width = 3, threshold = 0.05)
  expect_true(!is.null(result))
})

test_that("fregre.lime works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.lime(env$model, env$data, observation = 1,
                         n.samples = 50, kernel.width = 1, seed = 42)
  expect_true(!is.null(result))
})

test_that("fregre.counterfactual works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.counterfactual(env$model, env$data, observation = 1,
                                   target.value = 0.5)
  expect_true(!is.null(result))
})

test_that("fregre.anchor works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.anchor(env$model, env$data, observation = 1,
                           precision = 0.9, n.bins = 5)
  expect_true(!is.null(result))
})

test_that("fregre.friedman works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.friedman(env$model, env$data,
                             component.j = 1, component.k = 2, n.grid = 10)
  expect_true(!is.null(result))
})

test_that("fregre.depth works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.depth(env$model, env$data, env$y, n.boot = 10, seed = 42)
  expect_true(!is.null(result))
})

test_that("fregre.prototype works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.prototype(env$model, ncomp = 3,
                              n.prototypes = 3, n.criticisms = 3)
  expect_true(!is.null(result))
  expect_true("prototype_indices" %in% names(result) || "criticism_indices" %in% names(result))
})

test_that("fregre.stability works with logistic data", {
  env <- .test_logistic_model()
  result <- fregre.stability(env$data, env$y, ncomp = 3, n.boot = 10, seed = 42)
  expect_true(!is.null(result))
})

test_that("fregre.conditional.importance works with logistic", {
  env <- .test_logistic_model()
  result <- fregre.conditional.importance(env$model, env$data, env$y,
                                           n.bins = 3, n.perm = 10, seed = 42)
  expect_true(!is.null(result))
})

# =============================================================================
# Logistic-specific diagnostics
# =============================================================================

test_that("fregre.calibration works", {
  env <- .test_logistic_model()
  result <- fregre.calibration(env$model, env$y, n.groups = 5)
  expect_true(!is.null(result))
})

test_that("fregre.ece works", {
  env <- .test_logistic_model()
  result <- fregre.ece(env$model, env$y, n.bins = 5)
  expect_true(!is.null(result))
  expect_true(is.numeric(result) || is.list(result))
})

test_that("logistic explain input validation works", {
  expect_error(fregre.beta.decomp("not_a_model"))
  expect_error(fregre.pointwise("not_a_model"))
})
