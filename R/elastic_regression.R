#' Elastic Regression
#'
#' Alignment-integrated regression and classification for functional data.

#' Elastic Scalar-on-Function Regression
#'
#' Fits a functional linear model with simultaneous elastic alignment of the
#' functional covariate. The alignment and regression are jointly optimized.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param y Response vector (numeric).
#' @param ncomp.beta Number of basis functions for the beta coefficient (default 10).
#' @param lambda Regularization parameter (default 0).
#' @param max.iter Maximum iterations for the alignment-regression loop (default 20).
#' @param tol Convergence tolerance (default 1e-4).
#'
#' @return An object of class 'elastic.regression' with components:
#'   \item{alpha}{Intercept}
#'   \item{beta}{Beta coefficient function}
#'   \item{fitted.values}{Fitted response values}
#'   \item{residuals}{Residuals}
#'   \item{sse}{Sum of squared errors}
#'   \item{r.squared}{R-squared}
#'   \item{gammas}{Estimated warping functions}
#'   \item{aligned.srsfs}{Aligned SRSF transforms}
#'   \item{n.iter}{Number of iterations to convergence}
#'   \item{fdataobj}{Original functional data}
#'   \item{y}{Response vector}
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' fit <- elastic.regression(fd, y)
#' fit
#' }
#'
#' @export
elastic.regression <- function(fdataobj, y, ncomp.beta = 10, lambda = 0,
                               max.iter = 20, tol = 1e-4) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  n <- nrow(fdataobj$data)
  if (length(y) != n) {
    stop("Length of y must equal number of curves in fdataobj")
  }

  result <- .Call("wrap__elastic_regression_rust",
                  fdataobj$data, as.numeric(y),
                  as.numeric(fdataobj$argvals),
                  as.integer(ncomp.beta), as.numeric(lambda),
                  as.integer(max.iter), as.numeric(tol))

  if (is.null(result)) {
    stop("elastic.regression failed: check data dimensions and parameters")
  }

  beta_fd <- fdata(matrix(result$beta, nrow = 1), argvals = fdataobj$argvals,
                   names = list(main = "beta(t)", xlab = fdataobj$names$xlab,
                                ylab = "beta"))

  structure(
    list(
      alpha = result$alpha,
      beta = beta_fd,
      fitted.values = result$fitted_values,
      residuals = result$residuals,
      sse = result$sse,
      r.squared = result$r_squared,
      gammas = result$gammas,
      aligned.srsfs = result$aligned_srsfs,
      n.iter = result$n_iter,
      fdataobj = fdataobj,
      y = y,
      call = match.call()
    ),
    class = "elastic.regression"
  )
}

#' Elastic Logistic Classification
#'
#' Fits a functional logistic regression with elastic alignment.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param y Binary response vector (0/1).
#' @param ncomp.beta Number of basis functions for beta (default 10).
#' @param lambda Regularization parameter (default 0).
#' @param max.iter Maximum iterations (default 20).
#' @param tol Convergence tolerance (default 1e-4).
#'
#' @return An object of class 'elastic.logistic' with components:
#'   \item{alpha}{Intercept}
#'   \item{beta}{Beta coefficient function}
#'   \item{probabilities}{Predicted probabilities}
#'   \item{predicted.classes}{Predicted class labels (0/1)}
#'   \item{accuracy}{Classification accuracy}
#'   \item{loss}{Final loss value}
#'   \item{gammas}{Estimated warping functions}
#'   \item{aligned.srsfs}{Aligned SRSF transforms}
#'   \item{n.iter}{Number of iterations}
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- sample(0:1, 50, replace = TRUE)
#' fit <- elastic.logistic(fd, y)
#' fit
#' }
#'
#' @export
elastic.logistic <- function(fdataobj, y, ncomp.beta = 10, lambda = 0,
                             max.iter = 20, tol = 1e-4) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  n <- nrow(fdataobj$data)
  if (length(y) != n) {
    stop("Length of y must equal number of curves in fdataobj")
  }

  result <- .Call("wrap__elastic_logistic_rust",
                  fdataobj$data, as.numeric(y),
                  as.numeric(fdataobj$argvals),
                  as.integer(ncomp.beta), as.numeric(lambda),
                  as.integer(max.iter), as.numeric(tol))

  if (is.null(result)) {
    stop("elastic.logistic failed: check data dimensions and parameters")
  }

  beta_fd <- fdata(matrix(result$beta, nrow = 1), argvals = fdataobj$argvals,
                   names = list(main = "beta(t)", xlab = fdataobj$names$xlab,
                                ylab = "beta"))

  structure(
    list(
      alpha = result$alpha,
      beta = beta_fd,
      probabilities = result$probabilities,
      predicted.classes = result$predicted_classes,
      accuracy = result$accuracy,
      loss = result$loss,
      gammas = result$gammas,
      aligned.srsfs = result$aligned_srsfs,
      n.iter = result$n_iter,
      fdataobj = fdataobj,
      y = y,
      call = match.call()
    ),
    class = "elastic.logistic"
  )
}

#' Elastic Principal Component Regression
#'
#' Fits scalar-on-function regression using elastic FPCA scores as predictors.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param y Response vector (numeric).
#' @param ncomp Number of principal components (default 3).
#' @param pca.method PCA decomposition method: "vertical" (amplitude),
#'   "horizontal" (phase), or "joint" (combined).
#' @param lambda Regularization parameter for alignment (default 0).
#' @param max.iter Maximum alignment iterations (default 20).
#' @param tol Convergence tolerance (default 1e-4).
#'
#' @return An object of class 'elastic.pcr' with components:
#'   \item{alpha}{Intercept}
#'   \item{coefficients}{PC regression coefficients}
#'   \item{fitted.values}{Fitted response values}
#'   \item{sse}{Sum of squared errors}
#'   \item{r.squared}{R-squared}
#'   \item{pca.method}{PCA method used}
#'   \item{karcher.mean}{Karcher mean curve}
#'   \item{vert.scores}{Vertical FPCA scores (if applicable)}
#'   \item{horiz.scores}{Horizontal FPCA scores (if applicable)}
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' fit <- elastic.pcr(fd, y, ncomp = 2)
#' fit
#' }
#'
#' @export
elastic.pcr <- function(fdataobj, y, ncomp = 3,
                        pca.method = c("vertical", "horizontal", "joint"),
                        lambda = 0, max.iter = 20, tol = 1e-4) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  pca.method <- match.arg(pca.method)

  n <- nrow(fdataobj$data)
  if (length(y) != n) {
    stop("Length of y must equal number of curves in fdataobj")
  }

  result <- .Call("wrap__elastic_pcr_rust",
                  fdataobj$data, as.numeric(y),
                  as.numeric(fdataobj$argvals),
                  as.integer(ncomp), pca.method,
                  as.numeric(lambda),
                  as.integer(max.iter), as.numeric(tol))

  if (is.null(result)) {
    stop("elastic.pcr failed: check data dimensions and parameters")
  }

  structure(
    list(
      alpha = result$alpha,
      coefficients = result$coefficients,
      fitted.values = result$fitted_values,
      sse = result$sse,
      r.squared = result$r_squared,
      pca.method = pca.method,
      ncomp = ncomp,
      karcher.mean = result$karcher_mean,
      karcher.mean.srsf = result$karcher_mean_srsf,
      karcher.gammas = result$karcher_gammas,
      karcher.aligned = result$karcher_aligned,
      vert.scores = result$vert_scores,
      horiz.scores = result$horiz_scores,
      fdataobj = fdataobj,
      y = y,
      call = match.call()
    ),
    class = "elastic.pcr"
  )
}

#' @export
print.elastic.regression <- function(x, ...) {
  cat("Elastic Scalar-on-Function Regression\n")
  cat("  Intercept:", format(x$alpha, digits = 4), "\n")
  cat("  R-squared:", format(x$r.squared, digits = 4), "\n")
  cat("  Iterations:", x$n.iter, "\n")
  invisible(x)
}

#' @export
print.elastic.logistic <- function(x, ...) {
  cat("Elastic Logistic Classification\n")
  cat("  Accuracy:", format(x$accuracy * 100, digits = 1), "%\n")
  cat("  Iterations:", x$n.iter, "\n")
  invisible(x)
}

#' @export
print.elastic.pcr <- function(x, ...) {
  cat("Elastic Principal Component Regression\n")
  cat("  PCA method:", x$pca.method, "\n")
  cat("  Components:", x$ncomp, "\n")
  cat("  R-squared:", format(x$r.squared, digits = 4), "\n")
  invisible(x)
}

#' @export
plot.elastic.regression <- function(x, type = c("fit", "beta", "residuals"), ...) {
  type <- match.arg(type)
  if (type == "fit") {
    df <- data.frame(Observed = x$y, Fitted = x$fitted.values)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Observed, y = .data$Fitted)) +
      ggplot2::geom_point(color = "#2E8B57", size = 2.5, alpha = 0.7) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::labs(title = "Elastic Regression: Observed vs Fitted",
                    subtitle = paste("R\u00b2 =", round(x$r.squared, 3)),
                    x = "Observed", y = "Fitted") +
      ggplot2::theme_minimal()
    print(p)
  } else if (type == "beta") {
    plot(x$beta, main = "Elastic Beta Coefficient", ...)
  } else if (type == "residuals") {
    df <- data.frame(Fitted = x$fitted.values, Residuals = x$residuals)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Fitted, y = .data$Residuals)) +
      ggplot2::geom_point(color = "#2E8B57", size = 2, alpha = 0.7) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(title = "Elastic Regression: Residuals",
                    x = "Fitted Values", y = "Residuals") +
      ggplot2::theme_minimal()
    print(p)
  }
  invisible(x)
}

#' @export
plot.elastic.logistic <- function(x, type = c("fit", "beta", "probabilities"), ...) {
  type <- match.arg(type)
  if (type == "fit" || type == "probabilities") {
    df <- data.frame(
      Observation = seq_along(x$probabilities),
      Probability = x$probabilities,
      Class = factor(x$y, labels = c("0", "1"))
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Observation, y = .data$Probability,
                                           color = .data$Class)) +
      ggplot2::geom_point(size = 2.5, alpha = 0.7) +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
      ggplot2::scale_color_manual(values = c("0" = "#4A90D9", "1" = "#D55E00")) +
      ggplot2::labs(title = "Elastic Logistic: Predicted Probabilities",
                    subtitle = paste("Accuracy:", round(x$accuracy * 100, 1), "%"),
                    x = "Observation", y = "P(Class = 1)", color = "True Class") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")
    print(p)
  } else if (type == "beta") {
    plot(x$beta, main = "Elastic Logistic Beta Coefficient", ...)
  }
  invisible(x)
}

#' @export
plot.elastic.pcr <- function(x, type = c("fit", "coefficients", "residuals"), ...) {
  type <- match.arg(type)
  if (type == "fit") {
    df <- data.frame(Observed = x$y, Fitted = x$fitted.values)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Observed, y = .data$Fitted)) +
      ggplot2::geom_point(color = "#7B2D8E", size = 2.5, alpha = 0.7) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::labs(title = "Elastic PCR: Observed vs Fitted",
                    subtitle = paste("R\u00b2 =", round(x$r.squared, 3)),
                    x = "Observed", y = "Fitted") +
      ggplot2::theme_minimal()
    print(p)
  } else if (type == "coefficients") {
    df <- data.frame(
      PC = factor(paste0("PC", seq_along(x$coefficients)),
                  levels = paste0("PC", seq_along(x$coefficients))),
      Coefficient = x$coefficients
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PC, y = .data$Coefficient)) +
      ggplot2::geom_col(fill = "#7B2D8E", alpha = 0.7) +
      ggplot2::labs(title = "Elastic PCR Coefficients", x = NULL, y = "Coefficient") +
      ggplot2::theme_minimal()
    print(p)
  } else if (type == "residuals") {
    resid <- x$y - x$fitted.values
    df <- data.frame(Fitted = x$fitted.values, Residuals = resid)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Fitted, y = .data$Residuals)) +
      ggplot2::geom_point(color = "#7B2D8E", size = 2, alpha = 0.7) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(title = "Elastic PCR: Residuals",
                    x = "Fitted Values", y = "Residuals") +
      ggplot2::theme_minimal()
    print(p)
  }
  invisible(x)
}
