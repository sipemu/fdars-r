# Tests for explainability functions

# Helper: create a fitted fregre.lm model for testing
.test_fregre_model <- function(seed = 42) {
  set.seed(seed)
  t <- seq(0, 1, length.out = 30)
  n <- 25
  X <- matrix(0, n, 30)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1, 1, 0.3)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(30, sd = 0.1)
    y[i] <- 2 * amp + rnorm(1, sd = 0.2)
  }
  fd <- fdata(X, argvals = t)
  model <- fregre.lm(fd, y, ncomp = 3)
  list(model = model, data = fd, y = y)
}

test_that("fregre.beta.decomp returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.beta.decomp(env$model)
  expect_true(!is.null(result))
  expect_true("components" %in% names(result) || "coefficients" %in% names(result))
})

test_that("fregre.pointwise returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.pointwise(env$model)
  expect_true(!is.null(result))
  expect_true("importance" %in% names(result) || "importance_normalized" %in% names(result))
})

test_that("fregre.shap returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.shap(env$model, env$data)
  expect_true(!is.null(result))
  expect_true("values" %in% names(result) || "base_value" %in% names(result))
})

test_that("fregre.influence returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.influence(env$model, env$data)
  expect_true(!is.null(result))
  expect_true("leverage" %in% names(result) || "cooks_distance" %in% names(result))
})

test_that("fregre.vif returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.vif(env$model, env$data)
  expect_true(!is.null(result))
  expect_true("vif" %in% names(result))
})

test_that("fregre.dfbetas returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.dfbetas(env$model, env$data)
  expect_true(!is.null(result))
  expect_true("dfbetas" %in% names(result) || "dffits" %in% names(result))
})

test_that("fregre.saliency returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.saliency(env$model, env$data)
  expect_true(!is.null(result))
  expect_true("saliency_map" %in% names(result) || "mean_absolute_saliency" %in% names(result))
})

test_that("fregre.pdp returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.pdp(env$model, env$data, component = 1, n.grid = 10)
  expect_true(!is.null(result))
  expect_true("grid_values" %in% names(result) || "pdp_curve" %in% names(result))
})

test_that("fregre.ale returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.ale(env$model, env$data, component = 1, n.bins = 10)
  expect_true(!is.null(result))
  expect_true("ale_values" %in% names(result) || "bin_midpoints" %in% names(result))
})

test_that("fregre.importance returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.importance(env$model, env$data, env$y, n.perm = 20, seed = 42)
  expect_true(!is.null(result))
  expect_true("importance" %in% names(result) || "baseline_metric" %in% names(result))
})

test_that("fregre.loo returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.loo(env$model, env$data, env$y)
  expect_true(!is.null(result))
  expect_true("press" %in% names(result) || "loo_r_squared" %in% names(result))
})

test_that("fregre.sobol returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.sobol(env$model, env$data, env$y)
  expect_true(!is.null(result))
  expect_true("first_order" %in% names(result) || "total_order" %in% names(result))
})

test_that("fregre.domain returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.domain(env$model, window.width = 3, threshold = 0.05)
  expect_true(!is.null(result))
  expect_true("pointwise_importance" %in% names(result) || "intervals" %in% names(result))
})

test_that("fregre.prediction.interval returns correct structure", {
  env <- .test_fregre_model()
  new_fd <- fdata(env$data$data[1:3, , drop = FALSE], argvals = env$data$argvals)
  result <- fregre.prediction.interval(env$model, env$data, new_fd, confidence = 0.95)
  expect_true(!is.null(result))
  expect_true("predictions" %in% names(result) || "lower" %in% names(result))
})

test_that("fregre.lime returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.lime(env$model, env$data, observation = 1,
                         n.samples = 50, kernel.width = 1, seed = 42)
  expect_true(!is.null(result))
})

test_that("fregre.counterfactual returns correct structure", {
  env <- .test_fregre_model()
  target <- mean(env$y) + 1
  result <- fregre.counterfactual(env$model, env$data, observation = 1,
                                   target.value = target)
  expect_true(!is.null(result))
})

test_that("fregre.anchor returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.anchor(env$model, env$data, observation = 1,
                           precision = 0.9, n.bins = 5)
  expect_true(!is.null(result))
})

test_that("fregre.friedman returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.friedman(env$model, env$data,
                             component.j = 1, component.k = 2, n.grid = 10)
  expect_true(!is.null(result))
})

test_that("fregre.depth returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.depth(env$model, env$data, env$y, n.boot = 20, seed = 42)
  expect_true(!is.null(result))
})

test_that("fregre.prototype returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.prototype(env$model, ncomp = 3,
                              n.prototypes = 3, n.criticisms = 3)
  expect_true(!is.null(result))
  expect_true("prototype_indices" %in% names(result) || "criticism_indices" %in% names(result))
})

test_that("fregre.stability returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.stability(env$data, env$y, ncomp = 3, n.boot = 20, seed = 42)
  expect_true(!is.null(result))
})

test_that("fregre.conformal returns correct structure", {
  env <- .test_fregre_model()
  test_fd <- fdata(env$data$data[1:5, , drop = FALSE], argvals = env$data$argvals)
  result <- fregre.conformal(env$model, env$data, env$y, test_fd,
                              cal.fraction = 0.3, alpha = 0.1, seed = 42)
  expect_true(!is.null(result))
})

test_that("fregre.conditional.importance returns correct structure", {
  env <- .test_fregre_model()
  result <- fregre.conditional.importance(env$model, env$data, env$y,
                                           n.bins = 3, n.perm = 20, seed = 42)
  expect_true(!is.null(result))
})

test_that("fregre.significant.regions works", {
  env <- .test_fregre_model()
  result <- fregre.significant.regions(env$model, alpha = 0.05)
  # May return NULL if no significant regions
  expect_true(is.null(result) || is.list(result))
})

test_that("explain input validation works", {
  expect_error(fregre.beta.decomp("not_a_model"))
  expect_error(fregre.pointwise("not_a_model"))
})

test_that("elastic.attribution returns correct structure", {
  set.seed(42)
  t <- seq(0, 1, length.out = 20)
  n <- 20
  X <- matrix(0, n, 20)
  y <- numeric(n)
  for (i in 1:n) {
    amp <- rnorm(1, 1, 0.3)
    X[i, ] <- amp * sin(2 * pi * t) + rnorm(20, sd = 0.1)
    y[i] <- 2 * amp + rnorm(1, sd = 0.2)
  }
  fd <- fdata(X, argvals = t)

  result <- elastic.attribution(fd, y, ncomp = 2, pca.method = "vertical",
                                 max.iter = 5, n.perm = 20, seed = 42)
  expect_true(!is.null(result))
})
