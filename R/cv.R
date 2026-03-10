#' Unified K-Fold Cross-Validation for Functional Data
#'
#' Performs k-fold cross-validation for any functional data model, producing
#' out-of-fold (OOF) predictions where each observation is predicted exactly
#' once (when it is in the test fold). Supports both regression and
#' classification, with optional stratified fold assignment. Repeated
#' cross-validation (\code{nrep > 1}) runs multiple random fold partitions
#' to assess prediction variability.
#'
#' @param fdataobj An \code{fdata} object containing functional predictors.
#' @param y Response vector: numeric for regression, factor (or coerced to
#'   factor) for classification.
#' @param fit.fn A function with signature \code{function(fdataobj, y, ...)}
#'   that returns a fitted model object.
#' @param predict.fn Optional prediction function with signature
#'   \code{function(model, newdata, ...)}. If \code{NULL} (default), uses
#'   \code{predict(model, newdata)}.
#' @param kfold Number of folds (default 10).
#' @param nrep Number of repetitions (default 1). When \code{nrep > 1}, the
#'   entire k-fold procedure is repeated with different random fold assignments,
#'   producing an \code{n x nrep} matrix of predictions and per-observation
#'   variability estimates.
#' @param type One of \code{"regression"} or \code{"classification"}.
#'   Auto-detected from \code{y} if it is a factor.
#' @param stratified Logical; whether to stratify fold assignments so each fold
#'   has a similar distribution of \code{y} (default \code{TRUE}).
#' @param seed Optional integer seed for reproducibility of fold assignments.
#' @param ... Additional arguments passed through to \code{fit.fn}.
#'
#' @return An object of class \code{"cv.fdata"} with components:
#' \describe{
#'   \item{oof.predictions}{Out-of-fold predictions (numeric vector for
#'     regression, factor for classification). When \code{nrep > 1}, aggregated
#'     across repetitions (mean for regression, majority vote for
#'     classification).}
#'   \item{oof.probabilities}{For classification only: matrix of class
#'     probabilities if \code{predict.fn} returns them; otherwise \code{NULL}.}
#'   \item{folds}{Integer vector of fold assignments (from rep 1 for backward
#'     compatibility).}
#'   \item{fold.models}{List of per-fold fitted model objects. When
#'     \code{nrep > 1}, a list of \code{nrep} lists, each of length
#'     \code{kfold}.}
#'   \item{metrics}{Named list of overall performance metrics (computed on
#'     aggregated predictions when \code{nrep > 1}).}
#'   \item{fold.metrics}{Data frame of per-fold metrics (from rep 1 for backward
#'     compatibility).}
#'   \item{call}{The matched call.}
#'   \item{type}{Character: \code{"regression"} or \code{"classification"}.}
#'   \item{kfold}{Number of folds used.}
#'   \item{y}{The response vector.}
#'   \item{nrep}{Number of repetitions (only when \code{nrep > 1}).}
#'   \item{oof.matrix}{Matrix (\code{n x nrep}) of per-repetition predictions
#'     (only when \code{nrep > 1}).}
#'   \item{oof.sd}{Numeric vector of per-observation prediction SD across
#'     repetitions (only when \code{nrep > 1}). For classification, the
#'     proportion of reps disagreeing with the majority vote.}
#'   \item{folds.matrix}{Integer matrix (\code{n x nrep}) of fold assignments
#'     per repetition (only when \code{nrep > 1}).}
#'   \item{rep.metrics}{Data frame with one row per repetition containing
#'     per-rep metrics (only when \code{nrep > 1}).}
#'   \item{metrics.summary}{List with mean and sd of each metric across
#'     repetitions (only when \code{nrep > 1}).}
#' }
#'
#' @details
#' The \code{fit.fn} can internally perform nested cross-validation for
#' hyperparameter tuning. For example, passing a wrapper around
#' \code{\link{fregre.pc.cv}} as \code{fit.fn} gives proper nested CV: outer
#' folds produce unbiased OOF predictions while inner CV selects optimal
#' parameters on training data only.
#'
#' Each fold is wrapped in \code{tryCatch}: if a fold fails, its predictions
#' are set to \code{NA} and a warning is issued, but the remaining folds
#' continue.
#'
#' When \code{nrep > 1}, each repetition uses a different random fold partition.
#' If \code{seed} is provided, deterministic per-rep seeds are derived from it
#' for full reproducibility.
#'
#' @export
#' @examples
#' # Simple regression CV with fixed hyperparameters
#' set.seed(1)
#' t <- seq(0, 1, length.out = 30)
#' X <- matrix(0, 40, 30)
#' for (i in 1:40) X[i, ] <- sin(2*pi*t * i/40) + rnorm(30, sd = 0.1)
#' y <- rowMeans(X) + rnorm(40, sd = 0.1)
#' fd <- fdata(X, argvals = t)
#'
#' result <- cv.fdata(fd, y,
#'   fit.fn = function(fd, y, ...) fregre.pc(fd, y, ncomp = 3),
#'   kfold = 5, seed = 42)
#' print(result)
cv.fdata <- function(fdataobj, y, fit.fn, predict.fn = NULL,
                     kfold = 10, nrep = 1,
                     type = c("regression", "classification"),
                     stratified = TRUE, seed = NULL, ...) {
  cl <- match.call()

  # Validate inputs
  if (!inherits(fdataobj, "fdata")) {
    stop("'fdataobj' must be an fdata object")
  }
  n <- nrow(fdataobj$data)
  if (length(y) != n) {
    stop("Length of 'y' (", length(y), ") must equal number of curves (",
         n, ")")
  }
  if (!is.function(fit.fn)) {
    stop("'fit.fn' must be a function")
  }
  if (!is.null(predict.fn) && !is.function(predict.fn)) {
    stop("'predict.fn' must be a function or NULL")
  }
  kfold <- as.integer(kfold)
  if (kfold < 2 || kfold > n) {
    stop("'kfold' must be between 2 and n (", n, ")")
  }
  nrep <- as.integer(nrep)
  if (nrep < 1) {
    stop("'nrep' must be >= 1")
  }

  # Determine type
  if (is.factor(y)) {
    type <- "classification"
  } else {
    type <- match.arg(type)
  }
  if (type == "classification" && !is.factor(y)) {
    y <- as.factor(y)
  }

  # Seed management for repeated CV
  if (!is.null(seed) && nrep > 1) {
    set.seed(seed)
    rep_seeds <- sample.int(.Machine$integer.max, size = nrep)
  } else {
    rep_seeds <- NULL
  }

  # --- Single-rep path (nrep == 1): unchanged behavior ---
  if (nrep == 1) {
    folds <- .create_folds(y, kfold, type, stratified, seed)
    cv_single <- .run_one_cv(fdataobj, y, fit.fn, predict.fn, folds,
                             kfold, type, ...)
    metrics <- .compute_cv_metrics(y, cv_single$oof, type)
    fold.metrics <- .compute_fold_metrics(y, cv_single$oof, folds, kfold, type)

    return(structure(
      list(
        oof.predictions = cv_single$oof,
        oof.probabilities = NULL,
        folds = folds,
        fold.models = cv_single$fold.models,
        metrics = metrics,
        fold.metrics = fold.metrics,
        call = cl,
        type = type,
        kfold = kfold,
        y = y
      ),
      class = "cv.fdata"
    ))
  }

  # --- Repeated CV path (nrep > 1) ---
  if (type == "regression") {
    oof_mat <- matrix(NA_real_, n, nrep)
  } else {
    oof_mat <- matrix(NA_character_, n, nrep)
  }
  folds_mat <- matrix(NA_integer_, n, nrep)
  all_fold_models <- vector("list", nrep)

  for (r in seq_len(nrep)) {
    rep_seed <- if (!is.null(rep_seeds)) rep_seeds[r] else NULL
    folds_r <- .create_folds(y, kfold, type, stratified, rep_seed)
    folds_mat[, r] <- folds_r

    cv_r <- .run_one_cv(fdataobj, y, fit.fn, predict.fn, folds_r,
                        kfold, type, ...)
    oof_mat[, r] <- cv_r$oof
    all_fold_models[[r]] <- cv_r$fold.models
  }

  # Aggregate predictions
  if (type == "regression") {
    oof_agg <- rowMeans(oof_mat, na.rm = TRUE)
    oof_sd <- apply(oof_mat, 1, sd, na.rm = TRUE)
  } else {
    # Majority vote for classification
    oof_agg_char <- apply(oof_mat, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) == 0) return(NA_character_)
      tab <- table(x)
      names(tab)[which.max(tab)]
    })
    oof_agg <- factor(oof_agg_char, levels = levels(y))
    # Disagreement proportion: fraction of reps differing from majority
    oof_sd <- apply(oof_mat, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) <= 1) return(0)
      tab <- table(x)
      1 - max(tab) / sum(tab)
    })
  }

  # Metrics on aggregated predictions
  metrics <- .compute_cv_metrics(y, oof_agg, type)
  # Fold metrics from rep 1 for backward compat
  fold.metrics <- .compute_fold_metrics(y, oof_mat[, 1], folds_mat[, 1],
                                        kfold, type)

  # Per-rep metrics
  rep_metrics <- .compute_rep_metrics(y, oof_mat, nrep, type)

  # Metrics summary
  metrics_summary <- .compute_metrics_summary(rep_metrics, type)

  structure(
    list(
      oof.predictions = oof_agg,
      oof.probabilities = NULL,
      folds = folds_mat[, 1],
      fold.models = all_fold_models,
      metrics = metrics,
      fold.metrics = fold.metrics,
      call = cl,
      type = type,
      kfold = kfold,
      y = y,
      nrep = nrep,
      oof.matrix = oof_mat,
      oof.sd = oof_sd,
      folds.matrix = folds_mat,
      rep.metrics = rep_metrics,
      metrics.summary = metrics_summary
    ),
    class = "cv.fdata"
  )
}

#' Run One CV Repetition
#'
#' @param fdataobj An fdata object.
#' @param y Response vector.
#' @param fit.fn Fitting function.
#' @param predict.fn Prediction function or NULL.
#' @param folds Integer vector of fold assignments.
#' @param kfold Number of folds.
#' @param type "regression" or "classification".
#' @param ... Additional arguments passed to fit.fn.
#' @return List with oof (predictions vector) and fold.models.
#' @keywords internal
.run_one_cv <- function(fdataobj, y, fit.fn, predict.fn, folds, kfold,
                        type, ...) {
  n <- length(y)
  if (type == "regression") {
    oof <- rep(NA_real_, n)
  } else {
    oof <- rep(NA_character_, n)
  }
  fold.models <- vector("list", kfold)
  fold.warnings <- character(0)

  for (k in seq_len(kfold)) {
    test_idx <- which(folds == k)
    train_idx <- which(folds != k)

    tryCatch({
      fit <- fit.fn(fdataobj[train_idx, ], y[train_idx], ...)
      fold.models[[k]] <- fit

      if (is.null(predict.fn)) {
        preds <- predict(fit, fdataobj[test_idx, ])
      } else {
        preds <- predict.fn(fit, fdataobj[test_idx, ])
      }

      oof[test_idx] <- if (type == "classification") {
        as.character(preds)
      } else {
        preds
      }
    }, error = function(e) {
      msg <- paste0("Fold ", k, " failed: ", conditionMessage(e))
      fold.warnings <<- c(fold.warnings, msg)
      warning(msg, call. = FALSE)
    })
  }

  # Convert classification predictions to factor with same levels
  if (type == "classification") {
    oof <- factor(oof, levels = levels(y))
  }

  if (length(fold.warnings) > 0) {
    warning(length(fold.warnings), " of ", kfold, " folds failed",
            call. = FALSE)
  }

  list(oof = oof, fold.models = fold.models)
}

#' Compute Per-Repetition Metrics
#' @keywords internal
.compute_rep_metrics <- function(y, oof_mat, nrep, type) {
  if (type == "regression") {
    rep_df <- data.frame(
      rep = seq_len(nrep),
      RMSE = numeric(nrep),
      MAE = numeric(nrep),
      R2 = numeric(nrep)
    )
    for (r in seq_len(nrep)) {
      valid <- !is.na(oof_mat[, r])
      if (any(valid)) {
        y_v <- y[valid]
        o_v <- oof_mat[valid, r]
        rep_df$RMSE[r] <- pred.RMSE(y_v, o_v)
        rep_df$MAE[r] <- pred.MAE(y_v, o_v)
        rep_df$R2[r] <- pred.R2(y_v, o_v)
      } else {
        rep_df$RMSE[r] <- NA_real_
        rep_df$MAE[r] <- NA_real_
        rep_df$R2[r] <- NA_real_
      }
    }
  } else {
    rep_df <- data.frame(
      rep = seq_len(nrep),
      accuracy = numeric(nrep)
    )
    for (r in seq_len(nrep)) {
      oof_r <- factor(oof_mat[, r], levels = levels(y))
      valid <- !is.na(oof_r)
      if (any(valid)) {
        rep_df$accuracy[r] <- mean(y[valid] == oof_r[valid])
      } else {
        rep_df$accuracy[r] <- NA_real_
      }
    }
  }
  rep_df
}

#' Compute Metrics Summary Across Repetitions
#' @keywords internal
.compute_metrics_summary <- function(rep_metrics, type) {
  if (type == "regression") {
    list(
      RMSE = c(mean = mean(rep_metrics$RMSE, na.rm = TRUE),
               sd = sd(rep_metrics$RMSE, na.rm = TRUE)),
      MAE = c(mean = mean(rep_metrics$MAE, na.rm = TRUE),
              sd = sd(rep_metrics$MAE, na.rm = TRUE)),
      R2 = c(mean = mean(rep_metrics$R2, na.rm = TRUE),
             sd = sd(rep_metrics$R2, na.rm = TRUE))
    )
  } else {
    list(
      accuracy = c(mean = mean(rep_metrics$accuracy, na.rm = TRUE),
                   sd = sd(rep_metrics$accuracy, na.rm = TRUE))
    )
  }
}

#' Create Stratified or Random Fold Assignments
#'
#' @param y Response vector.
#' @param kfold Number of folds.
#' @param type \code{"regression"} or \code{"classification"}.
#' @param stratified Logical.
#' @param seed Optional seed.
#' @return Integer vector of fold assignments.
#' @keywords internal
.create_folds <- function(y, kfold, type, stratified, seed) {
  n <- length(y)
  if (!is.null(seed)) set.seed(seed)

  if (!stratified) {
    return(sample(rep(seq_len(kfold), length.out = n)))
  }

  folds <- integer(n)

  if (type == "classification") {
    # Stratify by class level
    for (lev in levels(y)) {
      idx <- which(y == lev)
      folds[idx] <- sample(rep(seq_len(kfold), length.out = length(idx)))
    }
  } else {
    # Stratify by quantile bin for regression
    bins <- cut(y, breaks = kfold, labels = FALSE, include.lowest = TRUE)
    for (b in seq_len(kfold)) {
      idx <- which(bins == b)
      folds[idx] <- sample(rep(seq_len(kfold), length.out = length(idx)))
    }
  }

  folds
}

#' Compute Overall CV Metrics
#' @keywords internal
.compute_cv_metrics <- function(y, oof, type) {
  valid <- !is.na(oof)
  if (sum(valid) == 0) {
    if (type == "regression") {
      return(list(RMSE = NA_real_, MAE = NA_real_, R2 = NA_real_))
    } else {
      return(list(accuracy = NA_real_))
    }
  }

  if (type == "regression") {
    y_v <- y[valid]
    o_v <- oof[valid]
    list(
      RMSE = pred.RMSE(y_v, o_v),
      MAE = pred.MAE(y_v, o_v),
      R2 = pred.R2(y_v, o_v)
    )
  } else {
    y_v <- y[valid]
    o_v <- oof[valid]
    conf <- table(Predicted = o_v, Actual = y_v)
    acc <- sum(diag(conf)) / sum(conf)
    per_class <- diag(conf) / colSums(conf)
    per_class[is.nan(per_class)] <- NA_real_
    list(
      accuracy = acc,
      per.class.accuracy = per_class,
      confusion = conf
    )
  }
}

#' Compute Per-Fold CV Metrics
#' @keywords internal
.compute_fold_metrics <- function(y, oof, folds, kfold, type) {
  if (type == "regression") {
    fold_df <- data.frame(
      fold = seq_len(kfold),
      n = integer(kfold),
      RMSE = numeric(kfold),
      MAE = numeric(kfold),
      R2 = numeric(kfold)
    )
    for (k in seq_len(kfold)) {
      idx <- which(folds == k)
      fold_df$n[k] <- length(idx)
      valid <- !is.na(oof[idx])
      if (any(valid)) {
        y_k <- y[idx[valid]]
        o_k <- oof[idx[valid]]
        fold_df$RMSE[k] <- pred.RMSE(y_k, o_k)
        fold_df$MAE[k] <- pred.MAE(y_k, o_k)
        fold_df$R2[k] <- pred.R2(y_k, o_k)
      } else {
        fold_df$RMSE[k] <- NA_real_
        fold_df$MAE[k] <- NA_real_
        fold_df$R2[k] <- NA_real_
      }
    }
  } else {
    fold_df <- data.frame(
      fold = seq_len(kfold),
      n = integer(kfold),
      accuracy = numeric(kfold)
    )
    for (k in seq_len(kfold)) {
      idx <- which(folds == k)
      fold_df$n[k] <- length(idx)
      valid <- !is.na(oof[idx])
      if (any(valid)) {
        fold_df$accuracy[k] <- mean(y[idx[valid]] == oof[idx[valid]])
      } else {
        fold_df$accuracy[k] <- NA_real_
      }
    }
  }
  fold_df
}

#' Print Method for cv.fdata
#'
#' @param x A \code{cv.fdata} object.
#' @param ... Ignored.
#'
#' @export
print.cv.fdata <- function(x, ...) {
  has_rep <- !is.null(x$nrep) && x$nrep > 1

  if (has_rep) {
    cat("Repeated K-Fold Cross-Validation (cv.fdata)\n")
  } else {
    cat("K-Fold Cross-Validation (cv.fdata)\n")
  }
  cat("  Type:", x$type, "\n")
  if (has_rep) {
    total_fits <- x$kfold * x$nrep
    cat("  Folds:", x$kfold, "\u00d7", x$nrep, "repetitions",
        paste0("(", total_fits, " total fits)"), "\n")
  } else {
    cat("  Folds:", x$kfold, "\n")
  }
  n_valid <- sum(!is.na(x$oof.predictions))
  n_total <- length(x$oof.predictions)
  cat("  Observations:", n_total,
      if (n_valid < n_total) paste0("(", n_valid, " predicted)") else "",
      "\n")

  if (has_rep) {
    cat("\nAggregated metrics (on mean predictions):\n")
  } else {
    cat("\nOverall metrics:\n")
  }
  if (x$type == "regression") {
    cat("  RMSE:", format(x$metrics$RMSE, digits = 4), "\n")
    cat("  MAE: ", format(x$metrics$MAE, digits = 4), "\n")
    cat("  R2:  ", format(x$metrics$R2, digits = 4), "\n")

    if (has_rep) {
      ms <- x$metrics.summary
      cat("\nPer-repetition metrics (mean \u00b1 sd):\n")
      cat("  RMSE:", format(ms$RMSE["mean"], digits = 4), "\u00b1",
          format(ms$RMSE["sd"], digits = 4), "\n")
      cat("  MAE: ", format(ms$MAE["mean"], digits = 4), "\u00b1",
          format(ms$MAE["sd"], digits = 4), "\n")
      cat("  R2:  ", format(ms$R2["mean"], digits = 4), "\u00b1",
          format(ms$R2["sd"], digits = 4), "\n")
      cat("\nPrediction variability:\n")
      cat("  Mean per-observation SD:", format(mean(x$oof.sd, na.rm = TRUE),
          digits = 4), "\n")
      cat("  Max per-observation SD: ", format(max(x$oof.sd, na.rm = TRUE),
          digits = 4), "\n")
    } else {
      cat("\nPer-fold RMSE range: [",
          format(min(x$fold.metrics$RMSE, na.rm = TRUE), digits = 4), ", ",
          format(max(x$fold.metrics$RMSE, na.rm = TRUE), digits = 4), "]\n",
          sep = "")
    }
  } else {
    cat("  Accuracy:", format(x$metrics$accuracy, digits = 4), "\n")
    if (!is.null(x$metrics$per.class.accuracy)) {
      cat("  Per-class accuracy:\n")
      pca <- x$metrics$per.class.accuracy
      for (i in seq_along(pca)) {
        cat("    ", names(pca)[i], ":", format(pca[i], digits = 4), "\n")
      }
    }

    if (has_rep) {
      ms <- x$metrics.summary
      cat("\nPer-repetition accuracy (mean \u00b1 sd):\n")
      cat("  Accuracy:", format(ms$accuracy["mean"], digits = 4), "\u00b1",
          format(ms$accuracy["sd"], digits = 4), "\n")
      cat("\nPrediction variability:\n")
      cat("  Mean disagreement rate:", format(mean(x$oof.sd, na.rm = TRUE),
          digits = 4), "\n")
      cat("  Max disagreement rate: ", format(max(x$oof.sd, na.rm = TRUE),
          digits = 4), "\n")
    } else {
      cat("\nPer-fold accuracy range: [",
          format(min(x$fold.metrics$accuracy, na.rm = TRUE), digits = 4), ", ",
          format(max(x$fold.metrics$accuracy, na.rm = TRUE), digits = 4),
          "]\n", sep = "")
    }
  }

  invisible(x)
}

#' Plot Method for cv.fdata
#'
#' For regression: observed vs predicted scatter plot and per-fold residual
#' boxplots (side by side via patchwork if available, otherwise just the
#' scatter). For classification: confusion matrix heatmap.
#'
#' @param x A \code{cv.fdata} object.
#' @param ... Ignored.
#'
#' @return A \code{ggplot} object (or a patchwork composition).
#'
#' @export
plot.cv.fdata <- function(x, ...) {
  valid <- !is.na(x$oof.predictions)
  has_rep <- !is.null(x$nrep) && x$nrep > 1

  if (x$type == "regression") {
    if (has_rep) {
      # Repeated CV: pointrange + SD histogram
      df <- data.frame(
        Observed = x$y[valid],
        Predicted = as.numeric(x$oof.predictions[valid]),
        SD = x$oof.sd[valid]
      )

      p1 <- ggplot(df, aes(x = .data$Observed, y = .data$Predicted)) +
        geom_pointrange(
          aes(ymin = .data$Predicted - .data$SD,
              ymax = .data$Predicted + .data$SD),
          colour = "steelblue", alpha = 0.6, size = 0.3
        ) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                    colour = "red") +
        labs(title = paste0("OOF: Observed vs Predicted (",
                            x$nrep, " reps, \u00b11 SD)"),
             x = "Observed", y = "Predicted") +
        theme_minimal()

      df_sd <- data.frame(SD = x$oof.sd[valid])
      p2 <- ggplot(df_sd, aes(x = .data$SD)) +
        geom_histogram(fill = "steelblue", alpha = 0.7, bins = 20) +
        labs(title = "Prediction Variability",
             x = "Per-observation SD", y = "Count") +
        theme_minimal()

      if (requireNamespace("patchwork", quietly = TRUE)) {
        return(patchwork::wrap_plots(p1, p2))
      }
      return(p1)
    }

    # Single-rep regression: existing behavior
    df <- data.frame(
      Observed = x$y[valid],
      Predicted = x$oof.predictions[valid],
      Fold = factor(x$folds[valid]),
      Residual = x$y[valid] - x$oof.predictions[valid]
    )

    p1 <- ggplot(df, aes(x = .data$Observed, y = .data$Predicted)) +
      geom_point(colour = "steelblue", alpha = 0.6) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                  colour = "red") +
      labs(title = "OOF: Observed vs Predicted",
           x = "Observed", y = "Predicted") +
      theme_minimal()

    p2 <- ggplot(df, aes(x = .data$Fold, y = .data$Residual)) +
      geom_boxplot(fill = "lightblue", alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
      labs(title = "Per-Fold Residuals", x = "Fold", y = "Residual") +
      theme_minimal()

    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(patchwork::wrap_plots(p1, p2))
    }
    return(p1)
  }

  # Classification: confusion matrix heatmap (same for single and repeated)
  conf <- x$metrics$confusion
  if (is.null(conf)) {
    message("No confusion matrix available to plot")
    return(invisible(NULL))
  }

  conf_df <- as.data.frame(conf)
  names(conf_df) <- c("Predicted", "Actual", "Count")
  conf_df$Predicted <- factor(conf_df$Predicted)
  conf_df$Actual <- factor(conf_df$Actual)

  title_str <- if (has_rep) {
    paste0("Confusion Matrix (Accuracy: ",
           format(x$metrics$accuracy, digits = 3),
           ", ", x$nrep, " reps)")
  } else {
    paste0("Confusion Matrix (Accuracy: ",
           format(x$metrics$accuracy, digits = 3), ")")
  }

  p <- ggplot(conf_df, aes(x = .data$Actual, y = .data$Predicted,
                            fill = .data$Count)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = .data$Count), size = 5) +
    scale_fill_viridis_c(direction = -1, guide = "none") +
    labs(title = title_str, x = "Actual", y = "Predicted") +
    coord_equal() +
    theme_minimal()

  p
}
