#' Function-on-Scalar Regression
#'
#' Functions for regression where the response is a function and predictors
#' are scalar variables.

#' Function-on-Scalar Regression (Penalized)
#'
#' Fits a penalized function-on-scalar regression model where the response
#' is a function Y_i(t) and predictors are scalar: \code{Y_i(t) = mu(t) + sum_j x_ij beta_j(t) + eps_i(t)}.
#'
#' @param fdataobj An object of class 'fdata' (functional response).
#' @param predictors A matrix of scalar predictors (n x p).
#' @param lambda Smoothing/penalty parameter (default 0).
#'
#' @return An object of class 'fosr' with components:
#'   \item{intercept}{Intercept function mu(t) as fdata}
#'   \item{beta}{Coefficient functions beta_j(t) as fdata}
#'   \item{fitted}{Fitted functional values as fdata}
#'   \item{residuals}{Residual functions as fdata}
#'   \item{r.squared}{Global R-squared}
#'   \item{r.squared.t}{Pointwise R-squared}
#'   \item{gcv}{GCV criterion value}
#'   \item{lambda}{Smoothing parameter used}
#'
#' @seealso \code{\link{fosr.fpc}} for FPC-based FOSR,
#'   \code{\link{fanova}} for functional ANOVA
#'
#' @examples
#' \donttest{
#' # Functional response: 50 curves observed at 10 time points
#' Y <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' # Two scalar predictors
#' X <- cbind(rnorm(50), rnorm(50))
#' fit <- fosr(Y, predictors = X, lambda = 1)
#' fit
#' }
#'
#' @export
fosr <- function(fdataobj, predictors, lambda = 0) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  predictors <- as.matrix(predictors)
  n <- nrow(fdataobj$data)
  m <- ncol(fdataobj$data)

  if (nrow(predictors) != n) {
    stop("Number of rows in predictors must equal number of curves")
  }

  result <- .Call("wrap__fosr_rust", fdataobj$data, predictors, as.numeric(lambda))

  if (is.null(result)) {
    stop("fosr failed: check data dimensions")
  }

  intercept_fd <- fdata(matrix(result$intercept, nrow = 1), argvals = fdataobj$argvals,
                        names = list(main = "mu(t)", xlab = fdataobj$names$xlab,
                                     ylab = "mu"))

  beta_fd <- fdata(result$beta, argvals = fdataobj$argvals,
                   names = list(main = "beta(t)", xlab = fdataobj$names$xlab,
                                ylab = "beta"))

  fitted_fd <- fdata(result$fitted, argvals = fdataobj$argvals,
                     names = list(main = "Fitted", xlab = fdataobj$names$xlab,
                                  ylab = fdataobj$names$ylab))

  residuals_fd <- fdata(result$residuals, argvals = fdataobj$argvals,
                        names = list(main = "Residuals", xlab = fdataobj$names$xlab,
                                     ylab = fdataobj$names$ylab))

  structure(
    list(
      intercept = intercept_fd,
      beta = beta_fd,
      fitted = fitted_fd,
      residuals = residuals_fd,
      r.squared = result$r_squared,
      r.squared.t = result$r_squared_t,
      beta.se = if (!is.null(result$beta_se)) result$beta_se else NULL,
      lambda = result$lambda,
      gcv = result$gcv,
      fdataobj = fdataobj,
      predictors = predictors,
      call = match.call()
    ),
    class = "fosr"
  )
}

#' FPC-based Function-on-Scalar Regression
#'
#' Fits a function-on-scalar regression using functional principal components.
#'
#' @param fdataobj An object of class 'fdata' (functional response).
#' @param predictors A matrix of scalar predictors (n x p).
#' @param ncomp Number of FPC components (default 3).
#'
#' @return An object of class 'fosr'.
#'
#' @examples
#' \donttest{
#' Y <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' X <- cbind(rnorm(50), rnorm(50))
#' fit <- fosr.fpc(Y, predictors = X, ncomp = 2)
#' fit
#' }
#'
#' @export
fosr.fpc <- function(fdataobj, predictors, ncomp = 3) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  predictors <- as.matrix(predictors)
  n <- nrow(fdataobj$data)

  if (nrow(predictors) != n) {
    stop("Number of rows in predictors must equal number of curves")
  }

  result <- .Call("wrap__fosr_fpc_rust", fdataobj$data, predictors, as.integer(ncomp))

  if (is.null(result)) {
    stop("fosr.fpc failed: check data dimensions")
  }

  intercept_fd <- fdata(matrix(result$intercept, nrow = 1), argvals = fdataobj$argvals,
                        names = list(main = "mu(t)", xlab = fdataobj$names$xlab,
                                     ylab = "mu"))

  beta_fd <- fdata(result$beta, argvals = fdataobj$argvals,
                   names = list(main = "beta(t)", xlab = fdataobj$names$xlab,
                                ylab = "beta"))

  fitted_fd <- fdata(result$fitted, argvals = fdataobj$argvals,
                     names = list(main = "Fitted", xlab = fdataobj$names$xlab,
                                  ylab = fdataobj$names$ylab))

  residuals_fd <- fdata(result$residuals, argvals = fdataobj$argvals,
                        names = list(main = "Residuals", xlab = fdataobj$names$xlab,
                                     ylab = fdataobj$names$ylab))

  structure(
    list(
      intercept = intercept_fd,
      beta = beta_fd,
      fitted = fitted_fd,
      residuals = residuals_fd,
      r.squared = result$r_squared,
      r.squared.t = result$r_squared_t,
      ncomp = result$ncomp,
      fdataobj = fdataobj,
      predictors = predictors,
      call = match.call()
    ),
    class = "fosr"
  )
}

#' Functional ANOVA
#'
#' Tests for differences in mean functions across groups using a
#' permutation-based F-test.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param groups Integer vector of group labels.
#' @param n.perm Number of permutations (default 1000).
#'
#' @return An object of class 'fanova' with components:
#'   \item{group.means}{Group mean functions as fdata}
#'   \item{overall.mean}{Overall mean function}
#'   \item{f.statistic.t}{Pointwise F-statistic}
#'   \item{global.statistic}{Global test statistic}
#'   \item{p.value}{P-value from permutation test}
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' groups <- rep(1:2, each = 25)
#' res <- fanova(fd, groups = groups, n.perm = 100)
#' res
#' }
#'
#' @export
fanova <- function(fdataobj, groups, n.perm = 1000) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  n <- nrow(fdataobj$data)
  if (length(groups) != n) {
    stop("Length of groups must equal number of curves")
  }

  result <- .Call("wrap__fanova_rust", fdataobj$data,
                  as.integer(groups), as.integer(n.perm))

  if (is.null(result)) {
    stop("fanova failed: check data dimensions and groups")
  }

  group_means_fd <- fdata(result$group_means, argvals = fdataobj$argvals,
                          names = list(main = "Group Means",
                                       xlab = fdataobj$names$xlab,
                                       ylab = fdataobj$names$ylab))

  structure(
    list(
      group.means = group_means_fd,
      overall.mean = result$overall_mean,
      f.statistic.t = result$f_statistic_t,
      global.statistic = result$global_statistic,
      p.value = result$p_value,
      n.perm = result$n_perm,
      n.groups = result$n_groups,
      group.labels = result$group_labels,
      fdataobj = fdataobj,
      groups = groups,
      call = match.call()
    ),
    class = "fanova"
  )
}

#' Predict from Function-on-Scalar Regression
#'
#' @param object A fitted 'fosr' object.
#' @param new.predictors Matrix of new scalar predictors. If NULL, returns fitted.
#' @param ... Additional arguments (ignored).
#'
#' @return An fdata object of predicted functional values.
#'
#' @export
predict.fosr <- function(object, new.predictors = NULL, ...) {
  if (is.null(new.predictors)) {
    return(object$fitted)
  }
  new.predictors <- as.matrix(new.predictors)

  result <- .Call("wrap__predict_fosr_rust",
                  as.numeric(object$intercept$data[1, ]),
                  object$beta$data,
                  new.predictors,
                  if (!is.null(object$lambda)) object$lambda else 0)

  fdata(result, argvals = object$fdataobj$argvals,
        names = list(main = "Predicted", xlab = object$fdataobj$names$xlab,
                     ylab = object$fdataobj$names$ylab))
}

#' Print method for fosr objects
#' @param x A fosr object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object \code{x}.
#' @export
print.fosr <- function(x, ...) {
  cat("Function-on-Scalar Regression\n")
  cat("=============================\n")
  cat("  Number of observations:", nrow(x$fdataobj$data), "\n")
  cat("  Number of predictors:", ncol(x$predictors), "\n")
  cat("  Evaluation points:", ncol(x$fdataobj$data), "\n")
  cat("  R-squared:", round(x$r.squared, 4), "\n")
  if (!is.null(x$lambda)) {
    cat("  Lambda:", x$lambda, "\n")
  }
  if (!is.null(x$ncomp)) {
    cat("  FPC components:", x$ncomp, "\n")
  }
  invisible(x)
}

#' Print method for fanova objects
#' @param x A fanova object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object \code{x}.
#' @export
print.fanova <- function(x, ...) {
  cat("Functional ANOVA\n")
  cat("================\n")
  cat("  Number of groups:", x$n.groups, "\n")
  cat("  Number of observations:", nrow(x$fdataobj$data), "\n")
  cat("  Global F-statistic:", round(x$global.statistic, 4), "\n")
  cat("  P-value:", format.pval(x$p.value), "\n")
  cat("  Permutations:", x$n.perm, "\n")
  invisible(x)
}

#' Plot method for fosr objects
#' @param x A fosr object.
#' @param ... Additional arguments (ignored).
#' @return A ggplot object.
#' @export
plot.fosr <- function(x, ...) {
  beta_data <- x$beta$data
  p <- nrow(beta_data)
  m <- ncol(beta_data)
  argvals <- x$fdataobj$argvals

  df <- data.frame(
    t = rep(argvals, p),
    beta = as.vector(t(beta_data)),
    predictor = factor(rep(seq_len(p), each = m))
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$t, y = .data$beta,
                                    color = .data$predictor)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::labs(x = x$fdataobj$names$xlab %||% "t",
                  y = "beta(t)",
                  title = "Coefficient Functions",
                  color = "Predictor")
}

#' Plot method for fanova objects
#' @param x A fanova object.
#' @param ... Additional arguments (ignored).
#' @return A ggplot object.
#' @export
plot.fanova <- function(x, ...) {
  argvals <- x$fdataobj$argvals
  m <- length(argvals)
  ng <- x$n.groups

  means_data <- x$group.means$data
  df_means <- data.frame(
    t = rep(argvals, ng),
    value = as.vector(t(means_data)),
    group = factor(rep(seq_len(ng), each = m))
  )

  ggplot2::ggplot(df_means, ggplot2::aes(x = .data$t, y = .data$value,
                                          color = .data$group)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(x = x$fdataobj$names$xlab %||% "t",
                  y = x$fdataobj$names$ylab %||% "X(t)",
                  title = paste0("Group Means (p = ",
                                 format.pval(x$p.value), ")"),
                  color = "Group")
}

# =============================================================================
# 2D Function-on-Scalar Regression
# =============================================================================

#' 2D Function-on-Scalar Regression
#'
#' Fits a penalized regression model where the response is a 2D surface
#' Y_i(s,t) and predictors are scalar.
#'
#' @param fdataobj An fdata object where each row is a flattened 2D surface
#'   (m1 * m2 grid points).
#' @param predictors A matrix of scalar predictors (n x p).
#' @param argvals.s Grid points along the s dimension.
#' @param argvals.t Grid points along the t dimension.
#' @param lambda.s Smoothing penalty in s direction (default 0).
#' @param lambda.t Smoothing penalty in t direction (default 0).
#'
#' @return An object of class 'fosr.2d' with components:
#'   \item{intercept}{Intercept surface as fdata}
#'   \item{beta}{Coefficient surfaces as fdata}
#'   \item{fitted}{Fitted surfaces as fdata}
#'   \item{residuals}{Residual surfaces as fdata}
#'   \item{r.squared}{Global R-squared}
#'   \item{r.squared.pointwise}{Pointwise R-squared values}
#'   \item{lambda.s}{Penalty in s direction}
#'   \item{lambda.t}{Penalty in t direction}
#'   \item{gcv}{GCV criterion}
#'   \item{grid}{List with argvals.s, argvals.t, m1, m2}
#'
#' @examples
#' \donttest{
#' n <- 30; m1 <- 5; m2 <- 5
#' fd <- fdata(matrix(rnorm(n * m1 * m2), n, m1 * m2))
#' X <- matrix(rnorm(n * 2), n, 2)
#' fit <- fosr.2d(fd, X, seq(0, 1, length.out = m1), seq(0, 1, length.out = m2))
#' }
#'
#' @export
fosr.2d <- function(fdataobj, predictors, argvals.s, argvals.t,
                     lambda.s = 0, lambda.t = 0) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }
  predictors <- as.matrix(predictors)

  result <- .Call("wrap__fosr_2d_rust",
                  fdataobj$data, predictors,
                  as.numeric(argvals.s), as.numeric(argvals.t),
                  as.numeric(lambda.s), as.numeric(lambda.t))

  if (is.null(result)) {
    stop("fosr.2d failed: check data dimensions and grid compatibility")
  }

  m <- length(argvals.s) * length(argvals.t)
  argvals <- seq_len(m)

  structure(
    list(
      intercept = fdata(matrix(result$intercept, nrow = 1), argvals = argvals),
      beta = fdata(result$beta, argvals = argvals),
      fitted = fdata(result$fitted, argvals = argvals),
      residuals = fdata(result$residuals, argvals = argvals),
      r.squared = result$r_squared,
      r.squared.pointwise = result$r_squared_pointwise,
      beta.se = result$beta_se,
      lambda.s = result$lambda_s,
      lambda.t = result$lambda_t,
      gcv = result$gcv,
      grid = list(argvals.s = argvals.s, argvals.t = argvals.t,
                  m1 = result$grid_m1, m2 = result$grid_m2),
      fdataobj = fdataobj,
      predictors = predictors,
      call = match.call()
    ),
    class = "fosr.2d"
  )
}

#' @export
predict.fosr.2d <- function(object, new.predictors = NULL, ...) {
  if (is.null(new.predictors)) {
    return(object$fitted)
  }
  new.predictors <- as.matrix(new.predictors)

  result <- .Call("wrap__predict_fosr_2d_rust",
                  as.numeric(object$intercept$data[1, ]),
                  object$beta$data,
                  new.predictors,
                  object$grid$argvals.s,
                  object$grid$argvals.t,
                  object$lambda.s,
                  object$lambda.t)

  if (is.null(result)) {
    stop("predict.fosr.2d failed")
  }

  m <- object$grid$m1 * object$grid$m2
  fdata(result, argvals = seq_len(m))
}

#' @export
print.fosr.2d <- function(x, ...) {
  cat("2D Function-on-Scalar Regression\n")
  cat("  Grid:", x$grid$m1, "x", x$grid$m2, "\n")
  cat("  Predictors:", ncol(x$predictors), "\n")
  cat("  R-squared:", round(x$r.squared, 4), "\n")
  cat("  Lambda (s, t):", x$lambda.s, ",", x$lambda.t, "\n")
  invisible(x)
}
