#' Model Explainability
#'
#' Explainability tools for functional regression models (fregre.lm
#' and fregre.logistic). Includes global explanations (PDP, ALE, Sobol),
#' local explanations (SHAP, LIME, anchors, counterfactuals), diagnostics,
#' and calibration metrics for logistic models.

# Internal helper: extract model components for Rust explain wrappers
.fit_args <- function(model) {
  if (!inherits(model, "fregre.lm")) {
    stop("model must be of class 'fregre.lm'")
  }
  list(
    intercept = model$intercept,
    beta_t = as.numeric(model$beta.t$data),
    beta_se = if (!is.null(model$beta.se)) as.numeric(model$beta.se) else rep(0, ncol(model$fdataobj$data)),
    gamma = if (!is.null(model$gamma)) as.numeric(model$gamma) else numeric(0),
    fitted_values = as.numeric(model$fitted.values),
    residuals = as.numeric(model$residuals),
    r_squared = model$r.squared,
    r_squared_adj = model$r.squared.adj,
    std_errors = as.numeric(model$std.errors),
    ncomp = as.integer(model$ncomp),
    coefficients = as.numeric(model$coefficients),
    residual_se = model$residual.se,
    gcv = model$gcv,
    fpca_mean = as.numeric(model$.fpca_mean),
    fpca_rotation_data = as.numeric(model$.fpca_rotation),
    fpca_rotation_nrow = as.integer(nrow(model$.fpca_rotation)),
    fpca_rotation_ncol = as.integer(ncol(model$.fpca_rotation)),
    fpca_scores_data = as.numeric(model$.fpca_scores),
    fpca_scores_nrow = as.integer(nrow(model$.fpca_scores)),
    fpca_scores_ncol = as.integer(ncol(model$.fpca_scores))
  )
}

# Internal helper: build common args list for Rust explain wrappers
.explain_call <- function(fn_name, model, ...) {
  fa <- .fit_args(model)
  args <- c(
    list(fn_name),
    list(fa$intercept, fa$beta_t, fa$beta_se, fa$gamma,
         fa$fitted_values, fa$residuals,
         fa$r_squared, fa$r_squared_adj, fa$std_errors,
         fa$ncomp, fa$coefficients, fa$residual_se, fa$gcv,
         fa$fpca_mean, fa$fpca_rotation_data,
         fa$fpca_rotation_nrow, fa$fpca_rotation_ncol,
         fa$fpca_scores_data, fa$fpca_scores_nrow, fa$fpca_scores_ncol),
    list(...)
  )
  do.call(.Call, args)
}

# Internal helper: build common args list for logistic explain wrappers
.logistic_call <- function(fn_name, model, ...) {
  fa <- .logistic_args(model)
  args <- c(
    list(fn_name),
    list(fa$intercept, fa$beta_t, fa$beta_se, fa$gamma,
         fa$probabilities, fa$predicted_classes,
         fa$ncomp, fa$accuracy, fa$std_errors,
         fa$coefficients, fa$log_likelihood, fa$iterations,
         fa$fpca_mean, fa$fpca_rotation_data,
         fa$fpca_rotation_nrow, fa$fpca_rotation_ncol,
         fa$fpca_scores_data, fa$fpca_scores_nrow, fa$fpca_scores_ncol),
    list(...)
  )
  do.call(.Call, args)
}

#' Functional Partial Dependence Plot
#'
#' Computes the partial dependence of the response on a specific FPC score,
#' along with ICE (Individual Conditional Expectation) curves.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#' @param data An fdata object (the training data).
#' @param component Which FPC component (1-based).
#' @param n.grid Number of grid points for the PDP (default 20).
#'
#' @return A list with grid_values, pdp_curve, ice_curves, and component.
#' @export
fregre.pdp <- function(model, data, component, n.grid = 20) {
  if (inherits(model, "fregre.logistic")) {
    .logistic_call("wrap__explain_pdp_logistic_rust", model,
                   data$data, as.integer(component - 1), as.integer(n.grid))
  } else {
    .explain_call("wrap__explain_pdp_rust", model,
                  data$data, as.integer(component - 1), as.integer(n.grid))
  }
}

#' Beta Decomposition
#'
#' Decomposes the beta coefficient function into per-FPC contributions.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#'
#' @return A list with components (per-FPC beta contributions), coefficients,
#'   and variance_proportion.
#' @export
fregre.beta.decomp <- function(model) {
  if (inherits(model, "fregre.logistic")) {
    .logistic_call("wrap__explain_beta_decomposition_logistic_rust", model)
  } else {
    .explain_call("wrap__explain_beta_decomposition_rust", model)
  }
}

#' Pointwise Importance
#'
#' Computes importance of each domain point for prediction.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#'
#' @return A list with importance, importance_normalized,
#'   component_importance, and score_variance.
#' @export
fregre.pointwise <- function(model) {
  if (inherits(model, "fregre.logistic")) {
    .logistic_call("wrap__explain_pointwise_importance_logistic_rust", model)
  } else {
    .explain_call("wrap__explain_pointwise_importance_rust", model)
  }
}

#' Significant Regions
#'
#' Identifies domain regions where the beta coefficient is significantly
#' different from zero based on standard errors.
#'
#' @param model A fitted \code{fregre.lm} model.
#' @param alpha Significance level (default 0.05).
#'
#' @return A list with start_idx, end_idx, and direction vectors,
#'   or NULL if no significant regions found.
#' @export
fregre.significant.regions <- function(model, alpha = 0.05) {
  if (!inherits(model, "fregre.lm") && !inherits(model, "fregre.logistic")) {
    stop("model must be of class 'fregre.lm' or 'fregre.logistic'")
  }
  .Call("wrap__explain_significant_regions_rust",
        as.numeric(model$beta.t$data),
        as.numeric(model$beta.se),
        as.numeric(alpha))
}

#' FPC SHAP Values
#'
#' Computes Shapley values for the FPC scores.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#' @param data An fdata object (the training data).
#' @param n.samples Number of samples for logistic SHAP (default 100).
#' @param seed Random seed (used for logistic models).
#'
#' @return A list with values (SHAP matrix), base_value, and mean_scores.
#' @export
fregre.shap <- function(model, data, n.samples = 100, seed = NULL) {
  if (inherits(model, "fregre.logistic")) {
    if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
    .logistic_call("wrap__explain_shap_logistic_rust", model,
                   data$data, as.integer(n.samples), as.integer(seed))
  } else {
    .explain_call("wrap__explain_shap_rust", model, data$data)
  }
}

#' Influence Diagnostics
#'
#' Computes leverage, Cook's distance, and related diagnostics.
#'
#' @param model A fitted \code{fregre.lm} model.
#' @param data An fdata object (the training data).
#'
#' @return A list with leverage, cooks_distance, p, and mse.
#' @export
fregre.influence <- function(model, data) {
  .explain_call("wrap__explain_influence_rust", model, data$data)
}

#' Variance Inflation Factors
#'
#' Computes VIF for FPC scores to detect multicollinearity.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#' @param data An fdata object (the training data).
#'
#' @return A list with vif, labels, mean_vif, n_moderate, and n_severe.
#' @export
fregre.vif <- function(model, data) {
  if (inherits(model, "fregre.logistic")) {
    .logistic_call("wrap__explain_vif_logistic_rust", model, data$data)
  } else {
    .explain_call("wrap__explain_vif_rust", model, data$data)
  }
}

#' DFBETAS and DFFITS
#'
#' Computes DFBETAS (influence on each coefficient) and DFFITS
#' (influence on fitted values) for each observation.
#'
#' @param model A fitted \code{fregre.lm} model.
#' @param data An fdata object (the training data).
#'
#' @return A list with dfbetas (matrix), dffits, studentized_residuals,
#'   dfbetas_cutoff, and dffits_cutoff.
#' @export
fregre.dfbetas <- function(model, data) {
  .explain_call("wrap__explain_dfbetas_rust", model, data$data)
}

#' Functional Saliency Map
#'
#' Computes a saliency map showing which domain regions most influence
#' each observation's prediction.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#' @param data An fdata object (the training data, used for lm models).
#'
#' @return A list with saliency_map (matrix) and mean_absolute_saliency.
#' @export
fregre.saliency <- function(model, data = NULL) {
  if (inherits(model, "fregre.logistic")) {
    .logistic_call("wrap__explain_saliency_logistic_rust", model)
  } else {
    .explain_call("wrap__explain_saliency_rust", model, data$data)
  }
}

#' FPC Permutation Importance
#'
#' Measures importance of each FPC by permuting scores and measuring
#' the increase in prediction error.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#' @param data An fdata object (the training data).
#' @param y Response vector.
#' @param n.perm Number of permutations (default 100).
#' @param seed Random seed.
#'
#' @return A list with importance, baseline_metric, and permuted_metric.
#' @export
fregre.importance <- function(model, data, y, n.perm = 100, seed = NULL) {
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  if (inherits(model, "fregre.logistic")) {
    .logistic_call("wrap__explain_importance_logistic_rust", model,
                   data$data, as.numeric(y),
                   as.integer(n.perm), as.integer(seed))
  } else {
    .explain_call("wrap__explain_importance_rust", model,
                  data$data, as.numeric(y),
                  as.integer(n.perm), as.integer(seed))
  }
}

#' LOO-CV and PRESS
#'
#' Computes leave-one-out cross-validation statistics.
#'
#' @param model A fitted \code{fregre.lm} model.
#' @param data An fdata object (the training data).
#' @param y Response vector.
#'
#' @return A list with press, loo_r_squared, loo_residuals, leverage, and tss.
#' @export
fregre.loo <- function(model, data, y) {
  .explain_call("wrap__explain_loo_rust", model,
                data$data, as.numeric(y))
}

#' Sobol Indices
#'
#' Computes first-order and total-order Sobol sensitivity indices for
#' each FPC component.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#' @param data An fdata object (the training data).
#' @param y Response vector (used for lm models).
#' @param n.samples Number of samples for logistic Sobol (default 500).
#' @param seed Random seed (used for logistic models).
#'
#' @return A list with first_order, total_order, var_y, and component_variance.
#' @export
fregre.sobol <- function(model, data, y = NULL, n.samples = 500, seed = NULL) {
  if (inherits(model, "fregre.logistic")) {
    if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
    .logistic_call("wrap__explain_sobol_logistic_rust", model,
                   data$data, as.integer(n.samples), as.integer(seed))
  } else {
    .explain_call("wrap__explain_sobol_rust", model,
                  data$data, as.numeric(y))
  }
}

#' Accumulated Local Effects
#'
#' Computes ALE for a specific FPC component.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#' @param data An fdata object.
#' @param component Which FPC component (1-based).
#' @param n.bins Number of ALE bins (default 20).
#'
#' @return A list with bin_midpoints, ale_values, bin_edges,
#'   bin_counts, and component.
#' @export
fregre.ale <- function(model, data, component, n.bins = 20) {
  if (inherits(model, "fregre.logistic")) {
    .logistic_call("wrap__explain_ale_logistic_rust", model,
                   data$data, as.integer(component - 1), as.integer(n.bins))
  } else {
    .explain_call("wrap__explain_ale_rust", model,
                  data$data, as.integer(component - 1), as.integer(n.bins))
  }
}

#' Domain Selection
#'
#' Identifies the most important domain intervals for prediction
#' using a sliding window on pointwise importance.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#' @param window.width Width of the sliding window (default 5).
#' @param threshold Minimum importance threshold (default 0.1).
#'
#' @return A list with pointwise_importance, intervals, window_width,
#'   and threshold.
#' @export
fregre.domain <- function(model, window.width = 5, threshold = 0.1) {
  if (inherits(model, "fregre.logistic")) {
    .logistic_call("wrap__explain_domain_logistic_rust", model,
                   as.integer(window.width), as.numeric(threshold))
  } else {
    .explain_call("wrap__explain_domain_rust", model,
                  as.integer(window.width), as.numeric(threshold))
  }
}

#' Prediction Intervals
#'
#' Computes prediction intervals for new observations.
#'
#' @param model A fitted \code{fregre.lm} model.
#' @param train.data An fdata object (training data).
#' @param new.data An fdata object (new observations for prediction).
#' @param confidence Confidence level (default 0.95).
#'
#' @return A list with predictions, lower, upper, prediction_se,
#'   confidence_level, t_critical, and residual_se.
#' @export
fregre.prediction.interval <- function(model, train.data, new.data,
                                       confidence = 0.95) {
  .explain_call("wrap__explain_prediction_intervals_rust", model,
                train.data$data, new.data$data, as.numeric(confidence))
}

#' LIME Explanation
#'
#' Computes a Local Interpretable Model-agnostic Explanation for a
#' single observation.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#' @param data An fdata object.
#' @param observation Index of the observation to explain (1-based).
#' @param n.samples Number of neighborhood samples (default 500).
#' @param kernel.width Kernel width for local weighting (default 1).
#' @param seed Random seed.
#'
#' @return A list with observation, attributions, local_intercept,
#'   local_r_squared, and kernel_width.
#' @export
fregre.lime <- function(model, data, observation, n.samples = 500,
                        kernel.width = 1, seed = NULL) {
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  if (inherits(model, "fregre.logistic")) {
    .logistic_call("wrap__explain_lime_logistic_rust", model,
                   data$data, as.integer(observation),
                   as.integer(n.samples), as.numeric(kernel.width),
                   as.integer(seed))
  } else {
    .explain_call("wrap__explain_lime_rust", model,
                  data$data, as.integer(observation),
                  as.integer(n.samples), as.numeric(kernel.width),
                  as.integer(seed))
  }
}

#' Counterfactual Explanation
#'
#' Finds the minimal change in FPC scores needed to achieve a target
#' prediction value. For logistic models, uses gradient-based optimization
#' to find the class-flipping counterfactual.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#' @param data An fdata object.
#' @param observation Index of the observation (1-based).
#' @param target.value Target prediction value (for lm models).
#' @param max.iter Maximum iterations for logistic counterfactual (default 100).
#' @param step.size Step size for logistic counterfactual (default 0.1).
#'
#' @return A list with observation, original_scores, counterfactual_scores,
#'   delta_scores, delta_function, distance, original_prediction,
#'   counterfactual_prediction, and found.
#' @export
fregre.counterfactual <- function(model, data, observation, target.value = NULL,
                                  max.iter = 100, step.size = 0.1) {
  if (inherits(model, "fregre.logistic")) {
    .logistic_call("wrap__explain_counterfactual_logistic_rust", model,
                   data$data, as.integer(observation),
                   as.integer(max.iter), as.numeric(step.size))
  } else {
    .explain_call("wrap__explain_counterfactual_rust", model,
                  data$data, as.integer(observation), as.numeric(target.value))
  }
}

#' Anchor Explanation
#'
#' Finds a rule (set of FPC score conditions) that anchors the prediction
#' to a narrow range with high precision.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#' @param data An fdata object.
#' @param observation Index of the observation (1-based).
#' @param precision Minimum precision threshold (default 0.95).
#' @param n.bins Number of bins for discretizing scores (default 10).
#'
#' @return A list with conditions, precision, coverage, n_matching,
#'   observation, and predicted_value.
#' @export
fregre.anchor <- function(model, data, observation, precision = 0.95,
                          n.bins = 10) {
  if (inherits(model, "fregre.logistic")) {
    .logistic_call("wrap__explain_anchor_logistic_rust", model,
                   data$data, as.integer(observation),
                   as.numeric(precision), as.integer(n.bins))
  } else {
    .explain_call("wrap__explain_anchor_rust", model,
                  data$data, as.integer(observation),
                  as.numeric(precision), as.integer(n.bins))
  }
}

#' Conformal Prediction
#'
#' Computes distribution-free prediction intervals using conformal inference.
#'
#' @param model A fitted \code{fregre.lm} model.
#' @param train.data An fdata object (training data).
#' @param train.y Training response vector.
#' @param test.data An fdata object (test data).
#' @param cal.fraction Fraction of training data for calibration (default 0.2).
#' @param alpha Miscoverage level (default 0.1 for 90% intervals).
#' @param seed Random seed.
#'
#' @return A list with predictions, lower, upper, residual_quantile,
#'   and coverage.
#' @export
fregre.conformal <- function(model, train.data, train.y, test.data,
                             cal.fraction = 0.2, alpha = 0.1, seed = NULL) {
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  fa <- .fit_args(model)
  .Call("wrap__explain_conformal_rust",
        fa$intercept, fa$beta_t, fa$beta_se, fa$gamma,
        fa$fitted_values, fa$residuals,
        fa$r_squared, fa$r_squared_adj, fa$std_errors,
        fa$ncomp, fa$coefficients, fa$residual_se, fa$gcv,
        fa$fpca_mean, fa$fpca_rotation_data,
        fa$fpca_rotation_nrow, fa$fpca_rotation_ncol,
        fa$fpca_scores_data, fa$fpca_scores_nrow, fa$fpca_scores_ncol,
        train.data$data, as.numeric(train.y),
        test.data$data,
        as.numeric(cal.fraction), as.numeric(alpha), as.integer(seed))
}

#' Explanation Stability
#'
#' Assesses how stable the model explanation is under bootstrap resampling.
#'
#' @param data An fdata object.
#' @param y Response vector.
#' @param ncomp Number of components.
#' @param n.boot Number of bootstrap resamples (default 100).
#' @param seed Random seed.
#'
#' @return A list with beta_t_std, coefficient_std, metric_std,
#'   beta_t_cv, importance_stability, and n_boot_success.
#' @export
fregre.stability <- function(data, y, ncomp, n.boot = 100, seed = NULL) {
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  .Call("wrap__explain_stability_rust",
        data$data, as.numeric(y),
        as.integer(ncomp), as.integer(n.boot), as.integer(seed))
}

#' Friedman H-Statistic
#'
#' Measures the interaction strength between two FPC components
#' using the Friedman H-statistic.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#' @param data An fdata object.
#' @param component.j First component index (1-based).
#' @param component.k Second component index (1-based).
#' @param n.grid Grid size for PDP computation (default 20).
#'
#' @return A list with component_j, component_k, h_squared,
#'   grid_j, grid_k, and pdp_2d.
#' @export
fregre.friedman <- function(model, data, component.j, component.k,
                            n.grid = 20) {
  if (inherits(model, "fregre.logistic")) {
    .logistic_call("wrap__explain_friedman_logistic_rust", model,
                   data$data, as.integer(component.j - 1),
                   as.integer(component.k - 1), as.integer(n.grid))
  } else {
    .explain_call("wrap__explain_friedman_rust", model,
                  data$data, as.integer(component.j - 1),
                  as.integer(component.k - 1), as.integer(n.grid))
  }
}

#' Regression Depth
#'
#' Computes regression depth for the fitted model and individual
#' observation score depths.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#' @param data An fdata object.
#' @param y Response vector.
#' @param n.boot Number of bootstrap resamples (default 100).
#' @param seed Random seed.
#'
#' @return A list with beta_depth, score_depths, mean_score_depth,
#'   and n_boot_success.
#' @export
fregre.depth <- function(model, data, y, n.boot = 100, seed = NULL) {
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  if (inherits(model, "fregre.logistic")) {
    .logistic_call("wrap__explain_regression_depth_logistic_rust", model,
                   data$data, as.numeric(y),
                   as.integer(n.boot), as.integer(seed))
  } else {
    .explain_call("wrap__explain_regression_depth_rust", model,
                  data$data, as.numeric(y),
                  as.integer(n.boot), as.integer(seed))
  }
}

#' Prototype and Criticism
#'
#' Finds representative (prototype) and atypical (criticism) observations
#' using the Maximum Mean Discrepancy (MMD) framework on FPCA scores.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model
#'   (uses its FPCA scores).
#' @param ncomp Number of components to use.
#' @param n.prototypes Number of prototypes (default 5).
#' @param n.criticisms Number of criticisms (default 5).
#'
#' @return A list with prototypes (indices), prototype_witness,
#'   criticisms (indices), criticism_witness, and bandwidth.
#' @export
fregre.prototype <- function(model, ncomp, n.prototypes = 5,
                             n.criticisms = 5) {
  if (!inherits(model, "fregre.lm") && !inherits(model, "fregre.logistic")) {
    stop("model must be of class 'fregre.lm' or 'fregre.logistic'")
  }
  .Call("wrap__explain_prototype_rust",
        as.numeric(model$.fpca_scores),
        as.integer(nrow(model$.fpca_scores)),
        as.integer(ncol(model$.fpca_scores)),
        as.integer(ncomp),
        as.integer(n.prototypes),
        as.integer(n.criticisms))
}

#' Calibration Diagnostics (Logistic)
#'
#' Computes calibration diagnostics for a functional logistic model:
#' Brier score, log loss, and Hosmer-Lemeshow test.
#'
#' @param model A fitted functional logistic model.
#' @param y True binary labels (0/1).
#' @param n.groups Number of groups for reliability diagram (default 10).
#'
#' @return A list with brier_score, log_loss, hosmer_lemeshow_chi2,
#'   hosmer_lemeshow_df, n_groups, reliability_bins, and bin_counts.
#' @export
fregre.calibration <- function(model, y, n.groups = 10) {
  # Extract logistic model components
  fa <- .logistic_args(model)
  .Call("wrap__explain_calibration_rust",
        fa$intercept, fa$beta_t, fa$beta_se, fa$gamma,
        fa$probabilities, fa$predicted_classes,
        fa$ncomp, fa$accuracy, fa$std_errors,
        fa$coefficients, fa$log_likelihood, fa$iterations,
        fa$fpca_mean, fa$fpca_rotation_data,
        fa$fpca_rotation_nrow, fa$fpca_rotation_ncol,
        fa$fpca_scores_data, fa$fpca_scores_nrow, fa$fpca_scores_ncol,
        as.numeric(y), as.integer(n.groups))
}

#' Expected Calibration Error (Logistic)
#'
#' Computes ECE, MCE, and ACE for a functional logistic model.
#'
#' @param model A fitted functional logistic model.
#' @param y True binary labels (0/1).
#' @param n.bins Number of bins (default 10).
#'
#' @return A list with ece, mce, ace, n_bins, and bin_ece_contributions.
#' @export
fregre.ece <- function(model, y, n.bins = 10) {
  fa <- .logistic_args(model)
  .Call("wrap__explain_ece_rust",
        fa$intercept, fa$beta_t, fa$beta_se, fa$gamma,
        fa$probabilities, fa$predicted_classes,
        fa$ncomp, fa$accuracy, fa$std_errors,
        fa$coefficients, fa$log_likelihood, fa$iterations,
        fa$fpca_mean, fa$fpca_rotation_data,
        fa$fpca_rotation_nrow, fa$fpca_rotation_ncol,
        fa$fpca_scores_data, fa$fpca_scores_nrow, fa$fpca_scores_ncol,
        as.numeric(y), as.integer(n.bins))
}

#' Conditional Permutation Importance
#'
#' Permutation importance respecting correlations between components.
#'
#' @param model A fitted \code{fregre.lm} or \code{fregre.logistic} model.
#' @param data An fdata object.
#' @param y Response vector.
#' @param n.bins Number of conditional bins (default 5).
#' @param n.perm Number of permutations (default 100).
#' @param seed Random seed.
#'
#' @return A list with importance, baseline_metric, and
#'   unconditional_importance.
#' @export
fregre.conditional.importance <- function(model, data, y, n.bins = 5,
                                          n.perm = 100, seed = NULL) {
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  if (inherits(model, "fregre.logistic")) {
    .logistic_call("wrap__explain_conditional_importance_logistic_rust", model,
                   data$data, as.numeric(y),
                   as.integer(n.bins), as.integer(n.perm), as.integer(seed))
  } else {
    .explain_call("wrap__explain_conditional_importance_rust", model,
                  data$data, as.numeric(y),
                  as.integer(n.bins), as.integer(n.perm), as.integer(seed))
  }
}

#' Elastic PCR Attribution
#'
#' Decomposes prediction importance into amplitude and phase contributions
#' after elastic PCR.
#'
#' @param fdataobj An fdata object.
#' @param y Response vector.
#' @param ncomp Number of elastic FPCA components.
#' @param pca.method PCA method: "vertical", "horizontal", or "joint".
#' @param lambda Regularization for alignment (default 0).
#' @param max.iter Maximum alignment iterations (default 20).
#' @param tol Convergence tolerance (default 1e-4).
#' @param n.perm Number of permutations (default 100).
#' @param seed Random seed.
#'
#' @return A list with amplitude_r_squared, phase_r_squared,
#'   total_r_squared, amplitude_importance, phase_importance,
#'   and p_values.
#' @export
elastic.attribution <- function(fdataobj, y, ncomp = 3,
                                pca.method = c("vertical", "horizontal", "joint"),
                                lambda = 0, max.iter = 20, tol = 1e-4,
                                n.perm = 100, seed = NULL) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  pca.method <- match.arg(pca.method)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  .Call("wrap__elastic_pcr_attribution_rust",
        fdataobj$data, as.numeric(y),
        as.numeric(fdataobj$argvals),
        as.integer(ncomp), pca.method,
        as.numeric(lambda), as.integer(max.iter), as.numeric(tol),
        as.integer(n.perm), as.integer(seed))
}

# Internal helper for logistic model components
.logistic_args <- function(model) {
  list(
    intercept = model$intercept,
    beta_t = as.numeric(model$beta.t$data),
    beta_se = if (!is.null(model$beta.se)) as.numeric(model$beta.se) else rep(0, ncol(model$fdataobj$data)),
    gamma = if (!is.null(model$gamma)) as.numeric(model$gamma) else numeric(0),
    probabilities = as.numeric(model$probabilities),
    predicted_classes = as.numeric(model$predicted.classes),
    ncomp = as.integer(model$ncomp),
    accuracy = model$accuracy,
    std_errors = if (!is.null(model$std.errors)) as.numeric(model$std.errors) else numeric(0),
    coefficients = as.numeric(model$coefficients),
    log_likelihood = model$log.likelihood,
    iterations = as.integer(model$iterations),
    fpca_mean = as.numeric(model$.fpca_mean),
    fpca_rotation_data = as.numeric(model$.fpca_rotation),
    fpca_rotation_nrow = as.integer(nrow(model$.fpca_rotation)),
    fpca_rotation_ncol = as.integer(ncol(model$.fpca_rotation)),
    fpca_scores_data = as.numeric(model$.fpca_scores),
    fpca_scores_nrow = as.integer(nrow(model$.fpca_scores)),
    fpca_scores_ncol = as.integer(ncol(model$.fpca_scores))
  )
}
