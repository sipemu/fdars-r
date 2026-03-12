#' Functional Mixed Models
#'
#' Functions for fitting functional mixed models with subject-level random effects.

#' Functional Mixed Model
#'
#' Fits a functional mixed model for repeated measures data:
#' \code{Y_ij(t) = mu(t) + sum_k x_ijk beta_k(t) + b_i(t) + eps_ij(t)}
#' where b_i(t) are subject-level random effects.
#'
#' @param fdataobj An object of class 'fdata'. All observations stacked
#'   (multiple per subject).
#' @param subject.ids Integer vector of subject identifiers (same length as nrow(data)).
#' @param covariates Optional matrix of covariates (n_total x p).
#' @param ncomp Number of FPC components for random effects (default 3).
#'
#' @return An object of class 'fmm' with components:
#'   \item{mean.function}{Overall mean function as fdata}
#'   \item{beta.functions}{Fixed effect coefficient functions as fdata}
#'   \item{random.effects}{Random effect functions per subject as fdata}
#'   \item{fitted}{Fitted values as fdata}
#'   \item{residuals}{Residuals as fdata}
#'   \item{sigma2.eps}{Residual variance estimate}
#'   \item{sigma2.u}{Random effect variance estimates}
#'   \item{n.subjects}{Number of unique subjects}
#'
#' @examples
#' \donttest{
#' # 10 subjects, 5 curves each = 50 total curves
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' subject <- rep(1:10, each = 5)
#' fit <- fmm(fd, subject.ids = subject)
#' fit
#' }
#'
#' @export
fmm <- function(fdataobj, subject.ids, covariates = NULL, ncomp = 3) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  n <- nrow(fdataobj$data)
  if (length(subject.ids) != n) {
    stop("Length of subject.ids must equal number of curves")
  }

  cov_mat <- if (!is.null(covariates)) as.matrix(covariates) else NULL

  result <- .Call("wrap__fmm_rust", fdataobj$data, as.integer(subject.ids),
                  cov_mat, as.integer(ncomp))

  if (is.null(result)) {
    stop("fmm failed: check data dimensions and subject IDs")
  }

  mean_fd <- fdata(matrix(result$mean_function, nrow = 1),
                   argvals = fdataobj$argvals,
                   names = list(main = "Mean Function",
                                xlab = fdataobj$names$xlab,
                                ylab = fdataobj$names$ylab))

  beta_fd <- fdata(result$beta_functions, argvals = fdataobj$argvals,
                   names = list(main = "Fixed Effects",
                                xlab = fdataobj$names$xlab,
                                ylab = "beta(t)"))

  random_fd <- fdata(result$random_effects, argvals = fdataobj$argvals,
                     names = list(main = "Random Effects",
                                  xlab = fdataobj$names$xlab,
                                  ylab = "b(t)"))

  fitted_fd <- fdata(result$fitted, argvals = fdataobj$argvals,
                     names = list(main = "Fitted", xlab = fdataobj$names$xlab,
                                  ylab = fdataobj$names$ylab))

  residuals_fd <- fdata(result$residuals, argvals = fdataobj$argvals,
                        names = list(main = "Residuals",
                                     xlab = fdataobj$names$xlab,
                                     ylab = fdataobj$names$ylab))

  structure(
    list(
      mean.function = mean_fd,
      beta.functions = beta_fd,
      random.effects = random_fd,
      fitted = fitted_fd,
      residuals = residuals_fd,
      random.variance = result$random_variance,
      sigma2.eps = result$sigma2_eps,
      sigma2.u = result$sigma2_u,
      ncomp = result$ncomp,
      n.subjects = result$n_subjects,
      eigenvalues = result$eigenvalues,
      fdataobj = fdataobj,
      subject.ids = subject.ids,
      covariates = covariates,
      call = match.call()
    ),
    class = "fmm"
  )
}

#' Predict from Functional Mixed Model
#'
#' Predict population-level (fixed effects only) functional responses
#' for new covariate values.
#'
#' @param object A fitted 'fmm' object.
#' @param new.covariates Matrix of new covariate values. If NULL, returns
#'   the overall mean function.
#' @param ... Additional arguments (ignored).
#'
#' @return An fdata object of predicted functional values.
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' subject <- rep(1:10, each = 5)
#' x <- cbind(rnorm(50))
#' fit <- fmm(fd, subject.ids = subject, covariates = x)
#' pred <- fmm.predict(fit, new.covariates = cbind(c(0, 1)))
#' }
#'
#' @export
fmm.predict <- function(object, new.covariates = NULL, ...) {
  cov_mat <- if (!is.null(new.covariates)) as.matrix(new.covariates) else NULL

  result <- .Call("wrap__fmm_predict_rust",
                  as.numeric(object$mean.function$data[1, ]),
                  object$beta.functions$data,
                  cov_mat)

  fdata(result, argvals = object$fdataobj$argvals,
        names = list(main = "Predicted",
                     xlab = object$fdataobj$names$xlab,
                     ylab = object$fdataobj$names$ylab))
}

#' Permutation Test for Fixed Effects in FMM
#'
#' Tests significance of fixed effects (covariates) using a
#' permutation-based F-test that respects the within-subject correlation.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param subject.ids Integer vector of subject identifiers.
#' @param covariates Matrix of covariates (n_total x p).
#' @param ncomp Number of FPC components (default 3).
#' @param n.perm Number of permutations (default 1000).
#' @param seed Random seed.
#'
#' @return An object of class 'fmm.test' with components:
#'   \item{f.statistics}{F-statistic per covariate}
#'   \item{p.values}{P-values per covariate}
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' subject <- rep(1:10, each = 5)
#' x <- cbind(rnorm(50))
#' test <- fmm.test.fixed(fd, subject.ids = subject, covariates = x,
#'                         n.perm = 100)
#' test
#' }
#'
#' @export
fmm.test.fixed <- function(fdataobj, subject.ids, covariates,
                             ncomp = 3, n.perm = 1000, seed = NULL) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  covariates <- as.matrix(covariates)
  seed_val <- if (!is.null(seed)) as.integer(seed) else 42L

  result <- .Call("wrap__fmm_test_fixed_rust", fdataobj$data,
                  as.integer(subject.ids), covariates,
                  as.integer(ncomp), as.integer(n.perm), seed_val)

  if (is.null(result)) {
    stop("fmm.test.fixed failed")
  }

  structure(
    list(
      f.statistics = result$f_statistics,
      p.values = result$p_values,
      n.perm = n.perm,
      call = match.call()
    ),
    class = "fmm.test"
  )
}

#' Print method for fmm objects
#' @param x A fmm object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object \code{x}.
#' @export
print.fmm <- function(x, ...) {
  cat("Functional Mixed Model\n")
  cat("======================\n")
  cat("  Number of observations:", nrow(x$fdataobj$data), "\n")
  cat("  Number of subjects:", x$n.subjects, "\n")
  cat("  FPC components:", x$ncomp, "\n")
  cat("  Residual variance (sigma2_eps):", round(x$sigma2.eps, 6), "\n")
  if (nrow(x$beta.functions$data) > 0) {
    cat("  Fixed effect covariates:", nrow(x$beta.functions$data), "\n")
  }
  invisible(x)
}

#' Print method for fmm.test objects
#' @param x A fmm.test object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object \code{x}.
#' @export
print.fmm.test <- function(x, ...) {
  cat("Fixed Effects Permutation Test\n")
  cat("==============================\n")
  cat("  Permutations:", x$n.perm, "\n\n")
  df <- data.frame(
    F.statistic = round(x$f.statistics, 4),
    P.value = format.pval(x$p.values)
  )
  rownames(df) <- paste0("Covariate ", seq_along(x$f.statistics))
  print(df)
  invisible(x)
}

#' Plot method for fmm objects
#' @param x A fmm object.
#' @param ... Additional arguments (ignored).
#' @return A ggplot object.
#' @export
plot.fmm <- function(x, ...) {
  argvals <- x$fdataobj$argvals
  m <- length(argvals)

  # Plot mean function and random effects
  n_subj <- x$n.subjects
  re_data <- x$random.effects$data

  df_re <- data.frame(
    t = rep(argvals, n_subj),
    value = as.vector(t(re_data)),
    subject = factor(rep(seq_len(n_subj), each = m))
  )

  df_mean <- data.frame(
    t = argvals,
    value = as.numeric(x$mean.function$data[1, ])
  )

  ggplot2::ggplot() +
    ggplot2::geom_line(data = df_re,
                       ggplot2::aes(x = .data$t, y = .data$value,
                                    group = .data$subject),
                       alpha = 0.3, color = "steelblue") +
    ggplot2::geom_line(data = df_mean,
                       ggplot2::aes(x = .data$t, y = .data$value),
                       linewidth = 1.2, color = "red") +
    ggplot2::labs(x = x$fdataobj$names$xlab %||% "t",
                  y = x$fdataobj$names$ylab %||% "X(t)",
                  title = "FMM: Mean (red) + Random Effects (blue)")
}
