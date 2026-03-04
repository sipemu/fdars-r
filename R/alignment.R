#' Elastic Alignment for Functional Data
#'
#' Functions for elastic curve alignment via SRSF transforms, Fisher-Rao distance,
#' and Karcher mean computation.

#' SRSF Transform
#'
#' Compute the Square-Root Slope Function (SRSF) transformation of functional data:
#' \eqn{q(t) = \text{sign}(f'(t)) \sqrt{|f'(t)|}}.
#'
#' @param fdataobj An object of class 'fdata'.
#'
#' @return An object of class 'fdata' containing the SRSF-transformed curves.
#'
#' @export
#' @examples
#' fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
#' q <- srsf.transform(fd)
srsf.transform <- function(fdataobj) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }
  result <- alignment_srsf_transform(fdataobj$data, fdataobj$argvals)
  fdata(result, argvals = fdataobj$argvals)
}

#' Inverse SRSF Transform
#'
#' Reconstruct a curve from its SRSF representation:
#' \eqn{f(t) = f_0 + \int_0^t q(s)|q(s)| ds}.
#'
#' @param q Numeric vector of SRSF values.
#' @param argvals Numeric vector of evaluation points.
#' @param f0 Initial value \eqn{f_0} (default 0).
#'
#' @return Numeric vector of reconstructed curve values.
#'
#' @export
#' @examples
#' fd <- fdata(matrix(sin(seq(0, pi, length.out = 50)), 1, 50),
#'             argvals = seq(0, 1, length.out = 50))
#' q <- srsf.transform(fd)
#' f_hat <- srsf.inverse(q$data[1, ], fd$argvals, fd$data[1, 1])
srsf.inverse <- function(q, argvals, f0 = 0) {
  alignment_srsf_inverse(as.numeric(q), as.numeric(argvals), as.numeric(f0))
}

#' Elastic Curve Alignment
#'
#' Align a set of functional curves using the elastic (Fisher-Rao) framework.
#' When no target is specified, the Karcher mean is used as the alignment target.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param target Optional target curve (numeric vector). If NULL, uses the
#'   cross-sectional mean as the alignment target.
#'
#' @return An object of class 'elastic.align' with components:
#' \describe{
#'   \item{aligned}{fdata of aligned curves}
#'   \item{gammas}{fdata of warping functions}
#'   \item{distances}{numeric vector of elastic distances}
#'   \item{target}{the target curve used for alignment}
#'   \item{fdataobj}{the original fdata input}
#' }
#'
#' @export
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
#' res <- elastic.align(fd)
#' }
elastic.align <- function(fdataobj, target = NULL) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  argvals <- fdataobj$argvals

  if (is.null(target)) {
    target <- colMeans(fdataobj$data)
  }

  res <- alignment_align_to_target(fdataobj$data, as.numeric(target), argvals)

  result <- list(
    aligned = fdata(res$aligned_data, argvals = argvals),
    gammas = fdata(res$gammas, argvals = argvals),
    distances = res$distances,
    target = target,
    fdataobj = fdataobj
  )
  class(result) <- "elastic.align"
  result
}

#' Elastic Distance Matrix
#'
#' Compute the pairwise elastic (Fisher-Rao) distance matrix for functional data.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param fdataref Optional reference 'fdata'. If NULL, computes self-distances.
#' @param ... Additional arguments (ignored).
#'
#' @return A distance matrix (numeric matrix).
#'
#' @export
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(100), 10, 10), argvals = seq(0, 1, length.out = 10))
#' D <- elastic.distance(fd)
#' }
elastic.distance <- function(fdataobj, fdataref = NULL, ...) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  if (is.null(fdataref)) {
    D <- alignment_self_dist(fdataobj$data, fdataobj$argvals)
  } else {
    if (!inherits(fdataref, "fdata")) {
      stop("fdataref must be of class 'fdata'")
    }
    D <- alignment_cross_dist(fdataobj$data, fdataref$data, fdataobj$argvals)
  }
  D
}

#' Elastic Distance (Metric Dispatcher Alias)
#'
#' Alias for \code{\link{elastic.distance}} for use with the \code{\link{metric}}
#' dispatcher.
#'
#' @inheritParams elastic.distance
#'
#' @return A distance matrix (numeric matrix).
#'
#' @export
metric.elastic <- function(fdataobj, fdataref = NULL, ...) {
  elastic.distance(fdataobj, fdataref, ...)
}

#' Karcher Mean in Elastic Metric
#'
#' Compute the Karcher (Fréchet) mean of functional data in the elastic metric.
#' This simultaneously estimates the mean shape and aligns all curves.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param max.iter Maximum number of iterations (default 20).
#' @param tol Convergence tolerance (default 1e-4).
#'
#' @return An object of class 'karcher.mean' with components:
#' \describe{
#'   \item{mean}{fdata of the Karcher mean curve}
#'   \item{aligned}{fdata of aligned curves}
#'   \item{gammas}{fdata of warping functions}
#'   \item{n.iter}{number of iterations used}
#'   \item{converged}{logical indicating convergence}
#'   \item{fdataobj}{the original fdata input}
#' }
#'
#' @export
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
#' km <- karcher.mean(fd)
#' }
karcher.mean <- function(fdataobj, max.iter = 20, tol = 1e-4) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  argvals <- fdataobj$argvals
  res <- alignment_karcher_mean(fdataobj$data, argvals,
                                as.integer(max.iter), as.numeric(tol))

  result <- list(
    mean = fdata(matrix(res$mean, nrow = 1), argvals = argvals),
    aligned = fdata(res$aligned_data, argvals = argvals),
    gammas = fdata(res$gammas, argvals = argvals),
    n.iter = res$n_iter,
    converged = res$converged,
    fdataobj = fdataobj
  )
  class(result) <- "karcher.mean"
  result
}

# =============================================================================
# S3 methods: elastic.align
# =============================================================================

#' @export
print.elastic.align <- function(x, ...) {
  n <- nrow(x$fdataobj$data)
  m <- ncol(x$fdataobj$data)
  cat("Elastic Alignment\n")
  cat("  Curves:", n, "x", m, "grid points\n")
  cat("  Mean elastic distance:", round(mean(x$distances), 4), "\n")
  invisible(x)
}

#' Plot Elastic Alignment Results
#'
#' @param x An object of class 'elastic.align'.
#' @param type Character: "aligned" (default), "original", "warps", or "both".
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object (invisibly).
#'
#' @export
plot.elastic.align <- function(x, type = c("aligned", "original", "warps", "both"), ...) {
  type <- match.arg(type)

  if (type == "both") {
    # Original and aligned side by side
    p1 <- .plot_fdata_curves(x$fdataobj, title = "Original")
    p2 <- .plot_fdata_curves(x$aligned, title = "Aligned")
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p <- p1 + p2
    } else {
      p <- p2
    }
  } else if (type == "warps") {
    p <- .plot_fdata_curves(x$gammas, title = "Warping Functions")
  } else if (type == "original") {
    p <- .plot_fdata_curves(x$fdataobj, title = "Original Curves")
  } else {
    p <- .plot_fdata_curves(x$aligned, title = "Aligned Curves")
  }

  print(p)
  invisible(p)
}

# =============================================================================
# S3 methods: karcher.mean
# =============================================================================

#' @export
print.karcher.mean <- function(x, ...) {
  n <- nrow(x$fdataobj$data)
  m <- ncol(x$fdataobj$data)
  cat("Karcher Mean (Elastic)\n")
  cat("  Curves:", n, "x", m, "grid points\n")
  cat("  Iterations:", x$n.iter, "\n")
  cat("  Converged:", x$converged, "\n")
  invisible(x)
}

#' Plot Karcher Mean Results
#'
#' @param x An object of class 'karcher.mean'.
#' @param type Character: "mean" (default), "aligned", or "warps".
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object (invisibly).
#'
#' @export
plot.karcher.mean <- function(x, type = c("mean", "aligned", "warps"), ...) {
  type <- match.arg(type)

  if (type == "aligned") {
    p <- .plot_fdata_curves(x$aligned, title = "Aligned Curves")
  } else if (type == "warps") {
    p <- .plot_fdata_curves(x$gammas, title = "Warping Functions")
  } else {
    # Mean curve with aligned curves in background
    n <- nrow(x$aligned$data)
    m <- ncol(x$aligned$data)
    argvals <- x$aligned$argvals

    df_aligned <- data.frame(
      curve_id = rep(seq_len(n), each = m),
      argval = rep(argvals, n),
      value = as.vector(t(x$aligned$data))
    )

    df_mean <- data.frame(
      argval = argvals,
      value = as.numeric(x$mean$data[1, ])
    )

    p <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = df_aligned,
        ggplot2::aes(x = .data$argval, y = .data$value, group = .data$curve_id),
        alpha = 0.3, color = "grey60"
      ) +
      ggplot2::geom_line(
        data = df_mean,
        ggplot2::aes(x = .data$argval, y = .data$value),
        color = "red", linewidth = 1.2
      ) +
      ggplot2::labs(title = "Karcher Mean", x = "t", y = "value") +
      ggplot2::theme_minimal()
  }

  print(p)
  invisible(p)
}

# =============================================================================
# Internal helper
# =============================================================================

#' @noRd
.plot_fdata_curves <- function(fd, title = "") {
  n <- nrow(fd$data)
  m <- ncol(fd$data)
  argvals <- fd$argvals

  df <- data.frame(
    curve_id = rep(seq_len(n), each = m),
    argval = rep(argvals, n),
    value = as.vector(t(fd$data))
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$argval, y = .data$value,
                                    group = .data$curve_id)) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::labs(title = title, x = "t", y = "value") +
    ggplot2::theme_minimal()
}
