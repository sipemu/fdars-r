#' Elastic Changepoint Detection
#'
#' Detect structural breaks in functional time series using elastic methods.

#' Elastic Changepoint Detection
#'
#' Tests for a changepoint in a functional time series using amplitude, phase,
#' or FPCA-based test statistics with permutation-based p-values.
#'
#' @param fdataobj An object of class 'fdata'. Curves should be in temporal
#'   order (rows = time points in the series).
#' @param type Type of changepoint test: "amplitude" (shape changes),
#'   "phase" (timing changes), or "fpca" (combined via elastic FPCA).
#' @param pca.method PCA method for type = "fpca": "vertical", "horizontal",
#'   or "joint".
#' @param ncomp Number of FPC components for type = "fpca" (default 3).
#' @param lambda Regularization for alignment (default 0).
#' @param max.iter Maximum alignment iterations (default 20).
#' @param n.mc Number of Monte Carlo permutations for p-value (default 1000).
#' @param cov.kernel Deprecated. No longer used (ignored with a warning).
#' @param cov.bandwidth Deprecated. No longer used (ignored with a warning).
#' @param seed Random seed for reproducibility.
#'
#' @return An object of class 'elastic.changepoint' with components:
#'   \item{changepoint}{Estimated changepoint index (1-based)}
#'   \item{test.statistic}{Value of the CUSUM test statistic}
#'   \item{p.value}{Permutation-based p-value}
#'   \item{cusum.values}{Full CUSUM process values}
#'   \item{type}{Type of test performed}
#'   \item{n}{Number of curves in the series}
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' cp <- elastic.changepoint(fd, type = "amplitude", n.mc = 100)
#' cp
#' }
#'
#' @export
elastic.changepoint <- function(fdataobj,
                                type = c("amplitude", "phase", "fpca"),
                                pca.method = c("vertical", "horizontal", "joint"),
                                ncomp = 3, lambda = 0, max.iter = 20,
                                n.mc = 1000,
                                cov.kernel = c("bartlett", "parzen", "flattop", "simple"),
                                cov.bandwidth = NULL,
                                seed = NULL) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  type <- match.arg(type)
  pca.method <- match.arg(pca.method)

  if (!missing(cov.kernel) && !is.null(cov.kernel)) {
    warning("'cov.kernel' is deprecated and ignored. Permutation-based p-values are used instead.")
  }
  if (!missing(cov.bandwidth) && !is.null(cov.bandwidth)) {
    warning("'cov.bandwidth' is deprecated and ignored. Permutation-based p-values are used instead.")
  }

  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  n <- nrow(fdataobj$data)
  argvals <- as.numeric(fdataobj$argvals)

  if (type == "amplitude") {
    result <- .Call("wrap__elastic_amp_changepoint_rust",
                    fdataobj$data, argvals,
                    as.numeric(lambda), as.integer(max.iter),
                    as.integer(n.mc), as.integer(seed))
  } else if (type == "phase") {
    result <- .Call("wrap__elastic_ph_changepoint_rust",
                    fdataobj$data, argvals,
                    as.numeric(lambda), as.integer(max.iter),
                    as.integer(n.mc), as.integer(seed))
  } else {
    result <- .Call("wrap__elastic_fpca_changepoint_rust",
                    fdataobj$data, argvals,
                    pca.method, as.integer(ncomp),
                    as.numeric(lambda), as.integer(max.iter),
                    as.integer(n.mc), as.integer(seed))
  }

  if (is.null(result)) {
    stop("elastic.changepoint failed: check data dimensions and parameters")
  }

  structure(
    list(
      changepoint = result$changepoint,
      test.statistic = result$test_statistic,
      p.value = result$p_value,
      cusum.values = result$cusum_values,
      type = type,
      n = n,
      fdataobj = fdataobj,
      call = match.call()
    ),
    class = "elastic.changepoint"
  )
}

#' @export
print.elastic.changepoint <- function(x, ...) {
  cat("Elastic Changepoint Detection\n")
  cat("  Type:", x$type, "\n")
  cat("  Changepoint at:", x$changepoint, "of", x$n, "curves\n")
  cat("  Test statistic:", format(x$test.statistic, digits = 4), "\n")
  cat("  p-value:", format(x$p.value, digits = 4), "\n")
  invisible(x)
}

#' @export
plot.elastic.changepoint <- function(x, type = c("statistic", "data"), ...) {
  type <- match.arg(type)
  if (type == "statistic") {
    df <- data.frame(
      Position = seq_along(x$cusum.values),
      CUSUM = x$cusum.values
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Position, y = .data$CUSUM)) +
      ggplot2::geom_line(color = "#333", linewidth = 1) +
      ggplot2::geom_vline(xintercept = x$changepoint, color = "#D55E00",
                          linewidth = 1, linetype = "dashed") +
      ggplot2::labs(title = paste0("Changepoint (", x$type, ")"),
                    subtitle = paste0("Detected at ", x$changepoint,
                                      " | p = ", format(x$p.value, digits = 3)),
                    x = "Position", y = "CUSUM Statistic") +
      ggplot2::theme_minimal()
    print(p)
  } else if (type == "data") {
    n <- x$n
    argvals <- x$fdataobj$argvals
    m <- length(argvals)
    df <- data.frame(
      t = rep(argvals, n),
      value = as.vector(t(x$fdataobj$data)),
      curve = factor(rep(seq_len(n), each = m)),
      segment = factor(rep(ifelse(seq_len(n) <= x$changepoint, "Before", "After"), each = m),
                       levels = c("Before", "After"))
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$t, y = .data$value,
                                           group = .data$curve, color = .data$segment)) +
      ggplot2::geom_line(alpha = 0.3) +
      ggplot2::scale_color_manual(values = c("Before" = "#4A90D9", "After" = "#D55E00")) +
      ggplot2::labs(title = paste0("Data with Changepoint at ", x$changepoint),
                    x = "t", y = "x(t)", color = "Segment") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")
    print(p)
  }
  invisible(x)
}
