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
#' @references
#' Srivastava, A., Klassen, E., Joshi, S.H., and Jermyn, I.H. (2011).
#' Shape analysis of elastic curves in Euclidean spaces.
#' \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence},
#' 33(7):1415--1428.
#'
#' Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
#' functional data using phase and amplitude separation.
#' \emph{Computational Statistics & Data Analysis}, 61:50--66.
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
#' @param periodic Logical; if TRUE, circularly rotate each curve to a canonical
#'   position before elastic alignment. This two-stage approach handles periodic
#'   functional data (e.g., data on \eqn{[0, 2\pi]} where \eqn{f(0) = f(2\pi)})
#'   that would otherwise be poorly aligned due to fixed boundary constraints
#'   \eqn{\gamma(0)=0, \gamma(1)=1}. Default is FALSE.
#' @param rotate.method Rotation method when \code{periodic = TRUE}: one of
#'   \code{"peak"} (default), \code{"xcorr"}, \code{"landmark"}, or
#'   \code{"iterative"}. See \code{\link{periodic.rotate}} for details.
#' @param rotate.args A named list of additional arguments passed to
#'   \code{\link{periodic.rotate}} (e.g., \code{reference}, \code{landmark.func},
#'   \code{max.iter}).
#'
#' @return An object of class 'elastic.align' with components:
#' \describe{
#'   \item{aligned}{fdata of aligned curves}
#'   \item{gammas}{fdata of warping functions}
#'   \item{distances}{numeric vector of elastic distances}
#'   \item{target}{the target curve used for alignment}
#'   \item{fdataobj}{the original fdata input}
#'   \item{rotations}{integer vector of circular rotation shifts applied
#'     (NULL when periodic = FALSE)}
#'   \item{rotate_method}{the rotation method used (NULL when periodic = FALSE)}
#' }
#'
#' @references
#' Srivastava, A., Klassen, E., Joshi, S.H., and Jermyn, I.H. (2011).
#' Shape analysis of elastic curves in Euclidean spaces.
#' \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence},
#' 33(7):1415--1428.
#'
#' Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
#' functional data using phase and amplitude separation.
#' \emph{Computational Statistics & Data Analysis}, 61:50--66.
#'
#' @export
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
#' res <- elastic.align(fd)
#' }
elastic.align <- function(fdataobj, target = NULL, periodic = FALSE,
                          rotate.method = "peak", rotate.args = list()) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  argvals <- fdataobj$argvals
  rotations <- NULL
  rotate_method_used <- NULL

  # Periodic pre-processing: circularly rotate curves to canonical position
  if (periodic) {
    rotate_method_used <- match.arg(rotate.method,
                                    c("peak", "xcorr", "landmark", "iterative"))

    # Extract reference from rotate.args if present
    ref_arg <- rotate.args$reference
    if (!is.null(ref_arg) && inherits(ref_arg, "fdata")) {
      ref_arg <- as.numeric(ref_arg$data[1, ])
    }

    rot <- do.call(.periodic_rotate, c(
      list(data_matrix = fdataobj$data,
           method = rotate_method_used,
           argvals = argvals),
      list(reference = ref_arg),
      rotate.args[!names(rotate.args) %in% "reference"]
    ))
    work_data <- rot$data
    rotations <- rot$shifts

    if (!is.null(target)) {
      if (inherits(target, "fdata")) {
        target <- as.numeric(target$data[1, ])
      }
      # Don't rotate the target when using xcorr/iterative — it IS the reference
      if (rotate_method_used %in% c("peak", "landmark")) {
        target_rot <- .periodic_rotate(matrix(target, nrow = 1),
                                       method = rotate_method_used,
                                       landmark_func = rotate.args$landmark.func)
        target <- as.numeric(target_rot$data[1, ])
      }
    } else {
      target <- colMeans(work_data)
    }
  } else {
    work_data <- fdataobj$data

    if (is.null(target)) {
      target <- colMeans(fdataobj$data)
    } else if (inherits(target, "fdata")) {
      target <- as.numeric(target$data[1, ])
    }
  }

  res <- alignment_align_to_target(work_data, as.numeric(target), argvals)

  result <- list(
    aligned = fdata(res$aligned_data, argvals = argvals),
    gammas = fdata(res$gammas, argvals = argvals),
    distances = res$distances,
    target = target,
    fdataobj = fdataobj,
    rotations = rotations,
    rotate_method = rotate_method_used
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
#' @references
#' Srivastava, A., Klassen, E., Joshi, S.H., and Jermyn, I.H. (2011).
#' Shape analysis of elastic curves in Euclidean spaces.
#' \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence},
#' 33(7):1415--1428.
#'
#' Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
#' functional data using phase and amplitude separation.
#' \emph{Computational Statistics & Data Analysis}, 61:50--66.
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
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
#' D <- metric.elastic(fd)
#' }
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
#' @param periodic Logical; if TRUE, circularly rotate each curve to a canonical
#'   position before computing the Karcher mean. See \code{\link{elastic.align}}
#'   for details. Default is FALSE.
#' @param rotate.method Rotation method when \code{periodic = TRUE}: one of
#'   \code{"peak"} (default), \code{"xcorr"}, \code{"landmark"}, or
#'   \code{"iterative"}. See \code{\link{periodic.rotate}} for details.
#' @param rotate.args A named list of additional arguments passed to
#'   \code{\link{periodic.rotate}} (e.g., \code{reference}, \code{landmark.func},
#'   \code{max.iter}).
#'
#' @return An object of class 'karcher.mean' with components:
#' \describe{
#'   \item{mean}{fdata of the Karcher mean curve}
#'   \item{aligned}{fdata of aligned curves}
#'   \item{gammas}{fdata of warping functions}
#'   \item{n.iter}{number of iterations used}
#'   \item{converged}{logical indicating convergence}
#'   \item{fdataobj}{the original fdata input}
#'   \item{rotations}{integer vector of circular rotation shifts applied
#'     (NULL when periodic = FALSE)}
#'   \item{rotate_method}{the rotation method used (NULL when periodic = FALSE)}
#' }
#'
#' @references
#' Srivastava, A., Klassen, E., Joshi, S.H., and Jermyn, I.H. (2011).
#' Shape analysis of elastic curves in Euclidean spaces.
#' \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence},
#' 33(7):1415--1428.
#'
#' Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
#' functional data using phase and amplitude separation.
#' \emph{Computational Statistics & Data Analysis}, 61:50--66.
#'
#' @export
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
#' km <- karcher.mean(fd)
#' }
karcher.mean <- function(fdataobj, max.iter = 20, tol = 1e-4, periodic = FALSE,
                         rotate.method = "peak", rotate.args = list()) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  argvals <- fdataobj$argvals
  rotations <- NULL
  rotate_method_used <- NULL

  if (periodic) {
    rotate_method_used <- match.arg(rotate.method,
                                    c("peak", "xcorr", "landmark", "iterative"))

    ref_arg <- rotate.args$reference
    if (!is.null(ref_arg) && inherits(ref_arg, "fdata")) {
      ref_arg <- as.numeric(ref_arg$data[1, ])
    }

    rot <- do.call(.periodic_rotate, c(
      list(data_matrix = fdataobj$data,
           method = rotate_method_used,
           argvals = argvals),
      list(reference = ref_arg),
      rotate.args[!names(rotate.args) %in% "reference"]
    ))
    work_data <- rot$data
    rotations <- rot$shifts
  } else {
    work_data <- fdataobj$data
  }

  res <- alignment_karcher_mean(work_data, argvals,
                                as.integer(max.iter), as.numeric(tol))

  result <- list(
    mean = fdata(matrix(res$mean, nrow = 1), argvals = argvals),
    mean_srsf = res$mean_srsf,
    aligned = fdata(res$aligned_data, argvals = argvals),
    gammas = fdata(res$gammas, argvals = argvals),
    n.iter = res$n_iter,
    converged = res$converged,
    fdataobj = fdataobj,
    rotations = rotations,
    rotate_method = rotate_method_used
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
  if (!is.null(x$rotations)) {
    cat("  Periodic: TRUE")
    if (!is.null(x$rotate_method)) {
      cat(" (method:", x$rotate_method, ")")
    }
    cat("\n")
  }
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

  p
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
  if (!is.null(x$rotations)) {
    cat("  Periodic: TRUE")
    if (!is.null(x$rotate_method)) {
      cat(" (method:", x$rotate_method, ")")
    }
    cat("\n")
  }
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

  p
}

# =============================================================================
# Periodic rotation
# =============================================================================

#' Periodic Rotation for Functional Data
#'
#' Circularly rotate each curve to a canonical position before elastic
#' alignment of periodic functional data (e.g., data on \eqn{[0, 2\pi]} where
#' \eqn{f(0) = f(2\pi)}).
#'
#' Four methods are available:
#' \describe{
#'   \item{\code{"peak"}}{Shift each curve's global maximum to \code{target.pos}
#'     using \code{which.max}. Fast and simple, but fragile when curves have
#'     multiple peaks of similar height.}
#'   \item{\code{"xcorr"}}{Circular cross-correlation with a reference curve
#'     (via \code{stats::convolve}). The lag maximising cross-correlation is the
#'     optimal shift. Robust to multi-peak shapes.}
#'   \item{\code{"landmark"}}{Like \code{"peak"}, but uses a user-supplied
#'     function \code{landmark.func(curve) -> index} to locate the feature of
#'     interest (e.g., steepest rise, first zero crossing).}
#'   \item{\code{"iterative"}}{Alternate cross-correlation rotation and elastic
#'     alignment for up to \code{max.iter} iterations. Useful when the optimal
#'     rotation depends on the alignment and vice versa.}
#' }
#'
#' @param fdataobj An object of class 'fdata'.
#' @param target.pos Target grid index for the feature position. If NULL
#'   (default), uses \code{floor(ncol(fdataobj$data) / 4)}.
#' @param method Rotation method: \code{"peak"} (default), \code{"xcorr"},
#'   \code{"landmark"}, or \code{"iterative"}.
#' @param reference Reference curve for cross-correlation methods. A numeric
#'   vector or single-curve \code{fdata}. If NULL, uses the cross-sectional
#'   mean.
#' @param landmark.func A function taking a numeric vector (one curve) and
#'   returning a single integer grid index. Required when
#'   \code{method = "landmark"}.
#' @param max.iter Maximum number of rotation-alignment iterations for
#'   \code{method = "iterative"}. Default 5.
#' @param ... Additional arguments (ignored).
#'
#' @return A list with components:
#' \describe{
#'   \item{fdataobj}{fdata of rotated curves}
#'   \item{shifts}{integer vector of circular shifts applied to each curve}
#' }
#'
#' @references
#' Srivastava, A. and Klassen, E. (2016). \emph{Functional and Shape Data
#' Analysis}. Springer.
#'
#' @export
#' @examples
#' # Periodic sinusoidal curves with random circular shifts
#' argvals <- seq(0, 2 * pi, length.out = 100)
#' data <- matrix(0, 10, 100)
#' for (i in 1:10) {
#'   shift <- sample(0:99, 1)
#'   idx <- ((seq_len(100) - 1 + shift) %% 100) + 1
#'   data[i, ] <- sin(argvals[idx])
#' }
#' fd <- fdata(data, argvals = argvals)
#' rot <- periodic.rotate(fd)
#'
#' # Cross-correlation method (robust to multi-peak curves)
#' rot_xcorr <- periodic.rotate(fd, method = "xcorr")
periodic.rotate <- function(fdataobj, target.pos = NULL,
                            method = c("peak", "xcorr", "landmark",
                                       "iterative"),
                            reference = NULL, landmark.func = NULL,
                            max.iter = 5, ...) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }
  method <- match.arg(method)

  if (!is.null(reference) && inherits(reference, "fdata")) {
    reference <- as.numeric(reference$data[1, ])
  }

  rot <- .periodic_rotate(fdataobj$data, target_pos = target.pos,
                          method = method, reference = reference,
                          landmark_func = landmark.func,
                          max_iter = max.iter,
                          argvals = fdataobj$argvals)
  list(
    fdataobj = fdata(rot$data, argvals = fdataobj$argvals),
    shifts = rot$shifts
  )
}

#' @noRd
.periodic_rotate <- function(data_matrix, target_pos = NULL,
                             method = "peak", reference = NULL,
                             landmark_func = NULL, max_iter = 5,
                             argvals = NULL) {
  method <- match.arg(method, c("peak", "xcorr", "landmark", "iterative"))

  switch(method,
    peak = .rotate_peak(data_matrix, target_pos),
    xcorr = .rotate_xcorr(data_matrix, reference),
    landmark = .rotate_landmark(data_matrix, target_pos, landmark_func),
    iterative = .rotate_iterative(data_matrix, max_iter, argvals)
  )
}

#' @noRd
.rotate_peak <- function(data_matrix, target_pos = NULL) {
  m <- ncol(data_matrix)
  n <- nrow(data_matrix)

  if (is.null(target_pos)) {
    target_pos <- floor(m / 4)
  }
  target_pos <- as.integer(target_pos)

  shifts <- integer(n)
  rotated <- matrix(0, n, m)

  for (i in seq_len(n)) {
    max_idx <- which.max(data_matrix[i, ])
    shift <- max_idx - target_pos
    shifts[i] <- as.integer(shift)
    idx <- ((seq_len(m) - 1L + shift) %% m) + 1L
    rotated[i, ] <- data_matrix[i, idx]
  }

  list(data = rotated, shifts = shifts)
}

#' @noRd
.rotate_xcorr <- function(data_matrix, reference = NULL) {
  m <- ncol(data_matrix)
  n <- nrow(data_matrix)

  if (is.null(reference)) {
    reference <- colMeans(data_matrix)
  }

  shifts <- integer(n)
  rotated <- matrix(0, n, m)

  for (i in seq_len(n)) {
    # Circular cross-correlation via FFT
    xcorr <- stats::convolve(data_matrix[i, ], reference, conj = TRUE,
                             type = "circular")
    # convolve with conj=TRUE computes sum_k x[k] * y[k+lag]
    # The maximum gives the best-matching lag
    best_lag <- which.max(xcorr) - 1L
    shifts[i] <- as.integer(best_lag)
    idx <- ((seq_len(m) - 1L + best_lag) %% m) + 1L
    rotated[i, ] <- data_matrix[i, idx]
  }

  list(data = rotated, shifts = shifts)
}

#' @noRd
.rotate_landmark <- function(data_matrix, target_pos = NULL,
                             landmark_func = NULL) {
  if (is.null(landmark_func)) {
    stop("landmark.func is required when method = 'landmark'")
  }

  m <- ncol(data_matrix)
  n <- nrow(data_matrix)

  if (is.null(target_pos)) {
    target_pos <- floor(m / 4)
  }
  target_pos <- as.integer(target_pos)

  shifts <- integer(n)
  rotated <- matrix(0, n, m)

  for (i in seq_len(n)) {
    lm_idx <- landmark_func(data_matrix[i, ])
    if (!is.numeric(lm_idx) || length(lm_idx) != 1 ||
        lm_idx < 1 || lm_idx > m || lm_idx != round(lm_idx)) {
      stop(sprintf(
        "landmark.func must return a single integer index in [1, %d]; got %s for curve %d",
        m, deparse(lm_idx), i
      ))
    }
    lm_idx <- as.integer(lm_idx)
    shift <- lm_idx - target_pos
    shifts[i] <- shift
    idx <- ((seq_len(m) - 1L + shift) %% m) + 1L
    rotated[i, ] <- data_matrix[i, idx]
  }

  list(data = rotated, shifts = shifts)
}

#' @noRd
.rotate_iterative <- function(data_matrix, max_iter = 5, argvals = NULL) {
  if (is.null(argvals)) {
    stop("argvals is required for method = 'iterative'")
  }

  n <- nrow(data_matrix)
  m <- ncol(data_matrix)
  current_data <- data_matrix
  cumulative_shifts <- integer(n)

  for (iter in seq_len(max_iter)) {
    # Step 1: compute reference as cross-sectional mean
    ref <- colMeans(current_data)

    # Step 2: cross-correlate all curves against reference
    rot <- .rotate_xcorr(current_data, reference = ref)

    # Accumulate shifts
    cumulative_shifts <- cumulative_shifts + rot$shifts

    # Step 3: check convergence (all shifts zero)
    if (all(rot$shifts == 0L)) break

    # Step 4: elastic alignment on rotated data
    aligned <- alignment_align_to_target(rot$data, ref, argvals)
    current_data <- aligned$aligned_data
  }

  # Apply cumulative shifts to original data for final output
  rotated <- matrix(0, n, m)
  for (i in seq_len(n)) {
    shift <- cumulative_shifts[i]
    idx <- ((seq_len(m) - 1L + shift) %% m) + 1L
    rotated[i, ] <- data_matrix[i, idx]
  }

  list(data = rotated, shifts = cumulative_shifts)
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

# =============================================================================
# Alignment Quality
# =============================================================================

#' Alignment Quality Diagnostics
#'
#' Compute comprehensive alignment quality metrics from a Karcher mean result,
#' including warp complexity, smoothness, and variance decomposition.
#'
#' @param fdataobj An object of class 'fdata' (the original data).
#' @param karcher An object of class 'karcher.mean'.
#'
#' @return An object of class 'alignment.quality' with components:
#' \describe{
#'   \item{warp_complexity}{Per-curve geodesic distance from warp to identity}
#'   \item{mean_warp_complexity}{Mean warp complexity}
#'   \item{warp_smoothness}{Per-curve bending energy}
#'   \item{mean_warp_smoothness}{Mean bending energy}
#'   \item{total_variance}{Total variance of original data}
#'   \item{amplitude_variance}{Amplitude variance (after alignment)}
#'   \item{phase_variance}{Phase variance (total - amplitude)}
#'   \item{phase_amplitude_ratio}{Phase-to-total variance ratio}
#'   \item{pointwise_variance_ratio}{Pointwise aligned/original variance ratio}
#'   \item{mean_variance_reduction}{Mean variance reduction}
#' }
#'
#' @export
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
#' km <- karcher.mean(fd, max.iter = 5)
#' aq <- alignment.quality(fd, km)
#' }
alignment.quality <- function(fdataobj, karcher) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }
  if (!inherits(karcher, "karcher.mean")) {
    stop("karcher must be of class 'karcher.mean'")
  }

  argvals <- fdataobj$argvals
  res <- alignment_quality_compute(
    fdataobj$data,
    as.numeric(karcher$mean$data[1, ]),
    as.numeric(alignment_srsf_transform(karcher$mean$data, argvals)[1, ]),
    karcher$gammas$data,
    karcher$aligned$data,
    argvals,
    as.integer(karcher$n.iter),
    karcher$converged
  )

  class(res) <- "alignment.quality"
  res
}

#' @export
print.alignment.quality <- function(x, ...) {
  cat("Alignment Quality Diagnostics\n")
  cat("  Mean warp complexity:", round(x$mean_warp_complexity, 4), "\n")
  cat("  Mean warp smoothness:", round(x$mean_warp_smoothness, 4), "\n")
  cat("  Total variance:     ", round(x$total_variance, 4), "\n")
  cat("  Amplitude variance: ", round(x$amplitude_variance, 4), "\n")
  cat("  Phase variance:     ", round(x$phase_variance, 4), "\n")
  cat("  Phase/Total ratio:  ", round(x$phase_amplitude_ratio, 4), "\n")
  cat("  Mean VR:            ", round(x$mean_variance_reduction, 4), "\n")
  invisible(x)
}

#' Plot Alignment Quality Diagnostics
#'
#' @param x An object of class 'alignment.quality'.
#' @param type Character: "variance" (default) or "complexity".
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object (invisibly).
#'
#' @export
plot.alignment.quality <- function(x, type = c("variance", "complexity"), ...) {
  type <- match.arg(type)

  if (type == "complexity") {
    df <- data.frame(
      curve = seq_along(x$warp_complexity),
      complexity = x$warp_complexity,
      smoothness = x$warp_smoothness
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$complexity, y = .data$smoothness)) +
      ggplot2::geom_point() +
      ggplot2::labs(title = "Warp Diagnostics",
                    x = "Warp Complexity", y = "Bending Energy") +
      ggplot2::theme_minimal()
  } else {
    df <- data.frame(
      component = c("Amplitude", "Phase"),
      variance = c(x$amplitude_variance, x$phase_variance)
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$component, y = .data$variance)) +
      ggplot2::geom_col(fill = c("steelblue", "coral")) +
      ggplot2::labs(title = "Variance Decomposition", x = "", y = "Variance") +
      ggplot2::theme_minimal()
  }

  p
}

# =============================================================================
# Elastic Decomposition
# =============================================================================

#' Elastic Phase-Amplitude Decomposition
#'
#' Decompose the difference between two curves into amplitude and phase
#' components using the elastic (Fisher-Rao) framework.
#'
#' @param f1 An 'fdata' object with a single curve, or a numeric vector.
#' @param f2 An 'fdata' object with a single curve, or a numeric vector.
#' @param argvals Evaluation points (required if f1/f2 are numeric vectors).
#' @param lambda Penalty weight on warp deviation from identity. Default 0.
#'
#' @return A list with components:
#' \describe{
#'   \item{gamma}{Warping function}
#'   \item{f_aligned}{Aligned version of f2}
#'   \item{d_amplitude}{Amplitude distance}
#'   \item{d_phase}{Phase distance}
#' }
#'
#' @export
#' @examples
#' \donttest{
#' t <- seq(0, 1, length.out = 50)
#' f1 <- fdata(matrix(sin(2*pi*t), 1, 50), argvals = t)
#' f2 <- fdata(matrix(sin(2*pi*t + 0.3), 1, 50), argvals = t)
#' dec <- elastic.decomposition(f1, f2)
#' }
elastic.decomposition <- function(f1, f2, argvals = NULL, lambda = 0) {
  if (inherits(f1, "fdata")) {
    argvals <- f1$argvals
    f1_vec <- as.numeric(f1$data[1, ])
  } else {
    f1_vec <- as.numeric(f1)
  }

  if (inherits(f2, "fdata")) {
    f2_vec <- as.numeric(f2$data[1, ])
  } else {
    f2_vec <- as.numeric(f2)
  }

  if (is.null(argvals)) {
    stop("argvals must be provided when f1/f2 are numeric vectors")
  }

  res <- alignment_decomposition(f1_vec, f2_vec, as.numeric(argvals),
                                 as.numeric(lambda))
  res
}

# =============================================================================
# Amplitude and Phase Distance
# =============================================================================

#' Amplitude Distance Matrix
#'
#' Compute the pairwise amplitude distance matrix. The amplitude distance
#' is the elastic (Fisher-Rao) distance after optimal alignment.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param fdataref Optional reference 'fdata'. If NULL, computes self-distances.
#' @param lambda Penalty weight. Default 0.
#' @param ... Additional arguments (ignored).
#'
#' @return A distance matrix.
#'
#' @export
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(60), 6, 10), argvals = seq(0, 1, length.out = 10))
#' D <- amplitude.distance(fd)
#' }
amplitude.distance <- function(fdataobj, fdataref = NULL, lambda = 0, ...) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  if (is.null(fdataref)) {
    D <- alignment_amplitude_dist(fdataobj$data, fdataobj$argvals,
                                  as.numeric(lambda))
  } else {
    if (!inherits(fdataref, "fdata")) {
      stop("fdataref must be of class 'fdata'")
    }
    # For cross-distances, use elastic cross-distance (amplitude = elastic distance)
    D <- alignment_cross_dist(fdataobj$data, fdataref$data, fdataobj$argvals)
  }

  D <- as.matrix(D)

  if (!is.null(rownames(fdataobj$data))) {
    rownames(D) <- rownames(fdataobj$data)
  }

  if (is.null(fdataref)) {
    if (!is.null(rownames(fdataobj$data))) {
      colnames(D) <- rownames(fdataobj$data)
    }
  } else if (!is.null(rownames(fdataref$data))) {
    colnames(D) <- rownames(fdataref$data)
  }

  D
}

#' Phase Distance Matrix
#'
#' Compute the pairwise phase distance matrix. The phase distance measures
#' the geodesic distance of the optimal warping function from the identity.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param fdataref Optional reference 'fdata'. If NULL, computes self-distances.
#' @param lambda Penalty weight. Default 0.
#' @param ... Additional arguments (ignored).
#'
#' @return A distance matrix.
#'
#' @export
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(60), 6, 10), argvals = seq(0, 1, length.out = 10))
#' D <- phase.distance(fd)
#' }
phase.distance <- function(fdataobj, fdataref = NULL, lambda = 0, ...) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  if (is.null(fdataref)) {
    D <- alignment_phase_dist(fdataobj$data, fdataobj$argvals,
                              as.numeric(lambda))
  } else {
    stop("phase.distance cross-distances not yet implemented; use self-distances")
  }

  D <- as.matrix(D)

  if (!is.null(rownames(fdataobj$data))) {
    rownames(D) <- rownames(fdataobj$data)
    colnames(D) <- rownames(fdataobj$data)
  }

  D
}

# =============================================================================
# Pairwise Consistency
# =============================================================================

#' Pairwise Alignment Consistency
#'
#' Measure the consistency of pairwise elastic alignments by checking the
#' triangle closure property. For each triplet (i, j, k), consistency checks
#' whether aligning i to k directly gives the same result as aligning i to j
#' then j to k.
#'
#' A value near 1 indicates high consistency (transitive alignments); a value
#' near 0 suggests the data contains distinct subgroups or the alignment is
#' sensitive to the target choice.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param lambda Penalty weight on warp deviation from identity. Default 0.
#' @param max.triplets Maximum number of random triplets to evaluate.
#'   Default 0 (use all triplets).
#'
#' @return A scalar in \[0, 1\] measuring alignment consistency.
#'
#' @export
#' @examples
#' \donttest{
#' t <- seq(0, 1, length.out = 100)
#' X <- matrix(0, 10, 100)
#' for (i in 1:10) X[i, ] <- sin(2*pi*(t - i/50))
#' fd <- fdata(X, argvals = t)
#' alignment.pairwise.consistency(fd, lambda = 0)
#' }
alignment.pairwise.consistency <- function(fdataobj, lambda = 0,
                                            max.triplets = 0) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be an fdata object")
  }

  alignment_pairwise_consistency(
    fdataobj$data,
    as.numeric(fdataobj$argvals),
    as.numeric(lambda),
    as.integer(max.triplets)
  )
}

# =============================================================================
# Constrained Alignment
# =============================================================================

#' Landmark-Constrained Elastic Alignment
#'
#' Align a set of curves to a target curve with landmark constraints.
#' The alignment is performed using segmented dynamic programming that
#' passes through specified landmark positions.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param target Target curve (numeric vector or single-curve fdata). If NULL,
#'   uses the cross-sectional mean.
#' @param landmark.pairs A two-column matrix of landmark pairs
#'   (target_t, source_t) for explicit constraints. If NULL, auto-detects.
#' @param kind Landmark type for auto-detection: "peak", "valley", "zero",
#'   "inflection". Default "peak".
#' @param min.prominence Minimum prominence for auto-detection. Default 0.
#' @param expected.count Expected number of landmarks for auto-detection.
#'   Default 0 (all detected).
#' @param lambda Penalty weight on warp deviation from identity. Default 0.
#'
#' @return An object of class 'elastic.align' with components:
#' \describe{
#'   \item{aligned}{fdata of aligned curves}
#'   \item{gammas}{fdata of warping functions}
#'   \item{distances}{numeric vector of elastic distances}
#'   \item{target}{the target curve}
#'   \item{fdataobj}{the original input}
#' }
#'
#' @export
#' @examples
#' \donttest{
#' t <- seq(0, 1, length.out = 100)
#' X <- matrix(0, 10, 100)
#' for (i in 1:10) X[i, ] <- sin(2*pi*(t - i/50))
#' fd <- fdata(X, argvals = t)
#' res <- elastic.align.constrained(fd, kind = "peak")
#' }
elastic.align.constrained <- function(fdataobj, target = NULL,
                                       landmark.pairs = NULL,
                                       kind = "peak", min.prominence = 0,
                                       expected.count = 0, lambda = 0) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  argvals <- fdataobj$argvals
  n <- nrow(fdataobj$data)

  if (is.null(target)) {
    target <- colMeans(fdataobj$data)
  } else if (inherits(target, "fdata")) {
    target <- as.numeric(target$data[1, ])
  }

  aligned_data <- matrix(0, n, length(argvals))
  gammas_data <- matrix(0, n, length(argvals))
  distances <- numeric(n)

  for (i in seq_len(n)) {
    fi <- fdataobj$data[i, ]

    if (!is.null(landmark.pairs)) {
      res <- alignment_constrained(
        as.numeric(target), as.numeric(fi), as.numeric(argvals),
        as.numeric(landmark.pairs[, 1]),
        as.numeric(landmark.pairs[, 2]),
        as.numeric(lambda)
      )
    } else {
      res <- alignment_with_landmarks(
        as.numeric(target), as.numeric(fi), as.numeric(argvals),
        kind, as.numeric(min.prominence),
        as.integer(expected.count), as.numeric(lambda)
      )
    }

    aligned_data[i, ] <- res$f_aligned
    gammas_data[i, ] <- res$gamma
    distances[i] <- res$distance
  }

  result <- list(
    aligned = fdata(aligned_data, argvals = argvals),
    gammas = fdata(gammas_data, argvals = argvals),
    distances = distances,
    target = target,
    fdataobj = fdataobj
  )
  class(result) <- "elastic.align"
  result
}

#' Apply Warping Function to a Curve
#'
#' Reparameterizes a curve f by composing it with a warping function gamma:
#' f(gamma(t)).
#'
#' @param fdataobj An fdata object (single curve, first row used).
#' @param gamma Numeric vector of warping function values on the same grid.
#'
#' @return An fdata object containing the reparameterized curve.
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(sin(seq(0, pi, length.out = 20)), 1, 20),
#'             argvals = seq(0, 1, length.out = 20))
#' gamma <- seq(0, 1, length.out = 20)^2
#' fd_warped <- srsf.reparameterize(fd, gamma)
#' }
#'
#' @export
srsf.reparameterize <- function(fdataobj, gamma) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }
  f <- as.numeric(fdataobj$data[1, ])
  argvals <- fdataobj$argvals
  gamma <- as.numeric(gamma)
  result <- alignment_reparameterize(f, argvals, gamma)
  fdata(matrix(result, nrow = 1), argvals = argvals, names = fdataobj$names)
}

#' Compose Two Warping Functions
#'
#' Computes the composition (gamma1 o gamma2)(t) = gamma1(gamma2(t)).
#'
#' @param gamma1 Numeric vector — outer warping function.
#' @param gamma2 Numeric vector — inner warping function.
#' @param fdataobj An fdata object providing the grid (argvals).
#'
#' @return Numeric vector of the composed warping function.
#'
#' @examples
#' \donttest{
#' t <- seq(0, 1, length.out = 20)
#' gamma1 <- t^2
#' gamma2 <- sqrt(t)
#' composed <- warp.compose(gamma1, gamma2, t)
#' }
#'
#' @export
warp.compose <- function(gamma1, gamma2, fdataobj) {
  if (inherits(fdataobj, "fdata")) {
    argvals <- fdataobj$argvals
  } else {
    argvals <- as.numeric(fdataobj)
  }
  as.numeric(alignment_compose_warps(as.numeric(gamma1),
                                     as.numeric(gamma2),
                                     argvals))
}

#' Elastic Pairwise Alignment
#'
#' Aligns curve f2 to curve f1 using elastic (Fisher-Rao) alignment.
#' Returns the aligned curve, warping function, and elastic distance.
#'
#' @param f1 An fdata object (reference curve, first row used).
#' @param f2 An fdata object (curve to align, first row used).
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{f.aligned} — fdata of the aligned curve.
#'   \item \code{gamma} — Numeric warping function.
#'   \item \code{distance} — Elastic distance after alignment.
#' }
#'
#' @examples
#' \donttest{
#' t <- seq(0, 1, length.out = 50)
#' f1 <- fdata(matrix(sin(2 * pi * t), 1, 50), argvals = t)
#' f2 <- fdata(matrix(sin(2 * pi * t^1.5), 1, 50), argvals = t)
#' res <- elastic.pair(f1, f2)
#' }
#'
#' @export
elastic.pair <- function(f1, f2) {
  if (!inherits(f1, "fdata") || !inherits(f2, "fdata")) {
    stop("f1 and f2 must be of class 'fdata'")
  }
  argvals <- f1$argvals
  res <- alignment_elastic_pair(as.numeric(f1$data[1, ]),
                                as.numeric(f2$data[1, ]),
                                as.numeric(argvals))
  list(
    f.aligned = fdata(matrix(res$f_aligned, nrow = 1), argvals = argvals,
                      names = f2$names),
    gamma = res$gamma,
    distance = res$distance
  )
}

#' Warping Complexity
#'
#' Computes the geodesic distance of a warping function from the identity,
#' measuring how much phase variability the warp introduces.
#'
#' @param gamma Numeric vector of warping function values.
#' @param argvals Numeric vector of grid points.
#'
#' @return A scalar complexity value (0 = identity warp).
#'
#' @examples
#' t <- seq(0, 1, length.out = 20)
#' warp.complexity(t^2, t)
#' warp.complexity(t, t)  # identity = 0
#'
#' @export
warp.complexity <- function(gamma, argvals) {
  alignment_warp_complexity(as.numeric(gamma), as.numeric(argvals))
}

#' Warping Smoothness
#'
#' Computes the bending energy of a warping function, measuring its roughness.
#'
#' @param gamma Numeric vector of warping function values.
#' @param argvals Numeric vector of grid points.
#'
#' @return A scalar smoothness (bending energy) value.
#'
#' @examples
#' t <- seq(0, 1, length.out = 20)
#' warp.smoothness(t^2, t)
#'
#' @export
warp.smoothness <- function(gamma, argvals) {
  alignment_warp_smoothness(as.numeric(gamma), as.numeric(argvals))
}

# =============================================================================
# Phase 6: Lambda CV, Warp Statistics, Phase Boxplots
# =============================================================================

#' Cross-Validation for Elastic Alignment Regularization Parameter
#'
#' Select the optimal regularization parameter (lambda) for elastic alignment
#' using K-fold cross-validation. The CV error measures how well the aligned
#' curves generalize to held-out folds.
#'
#' @param fdataobj An object of class \code{fdata}.
#' @param lambdas Numeric vector of candidate lambda values. Default is
#'   \code{10^seq(-4, 2, length.out = 20)}.
#' @param n.folds Number of cross-validation folds (default 5).
#' @param max.iter Maximum iterations for Karcher mean per fold (default 15).
#' @param tol Convergence tolerance (default 1e-3).
#' @param seed Random seed for fold assignment (default 42).
#'
#' @return An object of class \code{lambda.cv} with components:
#' \describe{
#'   \item{best.lambda}{The lambda value with the lowest CV error}
#'   \item{cv.scores}{Numeric vector of CV scores for each lambda}
#'   \item{lambdas}{The candidate lambda values tested}
#'   \item{call}{The matched call}
#' }
#'
#' @references
#' Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
#' functional data using phase and amplitude separation.
#' \emph{Computational Statistics & Data Analysis}, 61:50--66.
#'
#' @seealso \code{\link{karcher.mean}}, \code{\link{elastic.align}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 15, 50)
#' for (i in 1:15) X[i, ] <- sin(2 * pi * (t - i / 60))
#' fd <- fdata(X, argvals = t)
#' cv <- elastic.lambda.cv(fd, lambdas = 10^seq(-2, 1, length.out = 10))
#' cv
#' }
elastic.lambda.cv <- function(fdataobj,
                               lambdas = 10^seq(-4, 2, length.out = 20),
                               n.folds = 5, max.iter = 15, tol = 1e-3,
                               seed = 42) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  result <- .Call("wrap__alignment_lambda_cv_rust",
                  fdataobj$data, as.numeric(fdataobj$argvals),
                  as.numeric(lambdas), as.integer(n.folds),
                  as.integer(max.iter), as.numeric(tol), as.integer(seed))

  if (is.null(result)) {
    stop("elastic.lambda.cv failed: check data dimensions and parameters")
  }

  structure(
    list(
      best.lambda = result$best_lambda,
      cv.scores = result$cv_scores,
      lambdas = result$lambdas,
      call = match.call()
    ),
    class = "lambda.cv"
  )
}

#' Print Lambda CV Result
#'
#' @param x An object of class \code{lambda.cv}.
#' @param ... Additional arguments (ignored).
#' @export
print.lambda.cv <- function(x, ...) {
  cat("Lambda Cross-Validation\n")
  cat("  Candidates tested:", length(x$lambdas), "\n")
  cat("  Best lambda:", format(x$best.lambda, digits = 4), "\n")
  cat("  Best CV score:", format(min(x$cv.scores), digits = 4), "\n")
  invisible(x)
}

#' Plot Lambda CV Result
#'
#' Displays the cross-validation error as a function of the regularization
#' parameter lambda, with the optimal value marked.
#'
#' @param x An object of class \code{lambda.cv}.
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object.
#'
#' @export
plot.lambda.cv <- function(x, ...) {
  df <- data.frame(
    lambda = x$lambdas,
    cv_score = x$cv.scores
  )

  best_idx <- which.min(x$cv.scores)
  best_df <- data.frame(
    lambda = x$best.lambda,
    cv_score = x$cv.scores[best_idx]
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$lambda, y = .data$cv_score)) +
    ggplot2::geom_line(color = "#4A90D9", linewidth = 0.8) +
    ggplot2::geom_point(color = "#4A90D9", size = 1.5) +
    ggplot2::geom_point(
      data = best_df,
      ggplot2::aes(x = .data$lambda, y = .data$cv_score),
      color = "#D55E00", size = 3, shape = 17
    ) +
    ggplot2::geom_vline(
      xintercept = x$best.lambda,
      linetype = "dashed", color = "#D55E00", linewidth = 0.5
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = "Lambda Cross-Validation",
      x = expression(lambda), y = "CV Error"
    ) +
    ggplot2::theme_minimal()

  p
}

#' Warping Function Statistics
#'
#' Compute summary statistics for a collection of warping functions from
#' a Karcher mean computation, including the mean warp, variance, standard
#' deviation bands, and geodesic distances.
#'
#' @param karcher.result An object of class \code{karcher.mean}.
#'
#' @return An object of class \code{warp.statistics} with components:
#' \describe{
#'   \item{mean}{Mean warping function}
#'   \item{variance}{Pointwise variance of warps}
#'   \item{std.dev}{Pointwise standard deviation of warps}
#'   \item{lower.band}{Lower confidence band (mean - 2 SD)}
#'   \item{upper.band}{Upper confidence band (mean + 2 SD)}
#'   \item{karcher.mean.warp}{Karcher mean of warps in geodesic space}
#'   \item{geodesic.distances}{Geodesic distances of each warp from identity}
#'   \item{call}{The matched call}
#' }
#'
#' @references
#' Srivastava, A. and Klassen, E. (2016). \emph{Functional and Shape Data
#' Analysis}. Springer.
#'
#' @seealso \code{\link{karcher.mean}}, \code{\link{phase.boxplot}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 15, 50)
#' for (i in 1:15) X[i, ] <- sin(2 * pi * (t - i / 60))
#' fd <- fdata(X, argvals = t)
#' km <- karcher.mean(fd, max.iter = 5)
#' ws <- warp.statistics(km)
#' ws
#' }
warp.statistics <- function(karcher.result) {
  if (!inherits(karcher.result, "karcher.mean")) {
    stop("karcher.result must be of class 'karcher.mean'")
  }

  gammas <- karcher.result$gammas$data
  argvals <- as.numeric(karcher.result$gammas$argvals)

  result <- .Call("wrap__alignment_warp_statistics_rust",
                  gammas, argvals, 0.95)

  if (is.null(result)) {
    stop("warp.statistics failed: check karcher.result contents")
  }

  structure(
    list(
      mean = result$mean,
      variance = result$variance,
      std.dev = result$std_dev,
      lower.band = result$lower_band,
      upper.band = result$upper_band,
      karcher.mean.warp = result$karcher_mean_warp,
      geodesic.distances = result$geodesic_distances,
      argvals = argvals,
      call = match.call()
    ),
    class = "warp.statistics"
  )
}

#' Print Warp Statistics
#'
#' @param x An object of class \code{warp.statistics}.
#' @param ... Additional arguments (ignored).
#' @export
print.warp.statistics <- function(x, ...) {
  cat("Warping Function Statistics\n")
  cat("  Grid points:", length(x$mean), "\n")
  cat("  Curves:", length(x$geodesic.distances), "\n")
  cat("  Mean geodesic distance:", format(mean(x$geodesic.distances), digits = 4), "\n")
  cat("  Max geodesic distance:", format(max(x$geodesic.distances), digits = 4), "\n")
  invisible(x)
}

#' Plot Warp Statistics
#'
#' Displays the mean warping function with a shaded +/- 2 standard deviation
#' confidence band.
#'
#' @param x An object of class \code{warp.statistics}.
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object.
#'
#' @export
plot.warp.statistics <- function(x, ...) {
  argvals <- x$argvals
  df <- data.frame(
    t = argvals,
    mean = x$mean,
    lower = x$lower.band,
    upper = x$upper.band
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$t)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      fill = "#4A90D9", alpha = 0.3
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$mean),
      color = "#4A90D9", linewidth = 1
    ) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "grey50") +
    ggplot2::labs(
      title = "Warp Statistics",
      subtitle = "Mean warp \u00b1 2 SD",
      x = "t", y = expression(gamma(t))
    ) +
    ggplot2::theme_minimal()

  p
}

#' Phase Boxplot for Warping Functions
#'
#' Compute a functional boxplot for warping functions, identifying the
#' median, central 50% region, whiskers, and outliers based on functional
#' depth.
#'
#' @param karcher.result An object of class \code{karcher.mean}.
#' @param factor Whisker extension factor, analogous to the IQR multiplier
#'   in a classical boxplot (default 1.5).
#'
#' @return An object of class \code{phase.boxplot} with components:
#' \describe{
#'   \item{median}{Median warping function}
#'   \item{median.index}{Index of the median warp in the input}
#'   \item{central.lower}{Lower boundary of the central 50% region}
#'   \item{central.upper}{Upper boundary of the central 50% region}
#'   \item{whisker.lower}{Lower whisker}
#'   \item{whisker.upper}{Upper whisker}
#'   \item{depths}{Functional depth values for each warp}
#'   \item{outlier.indices}{Indices of outlying warps}
#'   \item{factor}{The whisker extension factor used}
#'   \item{argvals}{Grid points}
#'   \item{call}{The matched call}
#' }
#'
#' @references
#' Sun, Y. and Genton, M.G. (2011). Functional boxplots.
#' \emph{Journal of Computational and Graphical Statistics}, 20(2):316--334.
#'
#' @seealso \code{\link{karcher.mean}}, \code{\link{warp.statistics}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 15, 50)
#' for (i in 1:15) X[i, ] <- sin(2 * pi * (t - i / 60))
#' fd <- fdata(X, argvals = t)
#' km <- karcher.mean(fd, max.iter = 5)
#' pb <- phase.boxplot(km)
#' pb
#' }
phase.boxplot <- function(karcher.result, factor = 1.5) {
  if (!inherits(karcher.result, "karcher.mean")) {
    stop("karcher.result must be of class 'karcher.mean'")
  }

  gammas <- karcher.result$gammas$data
  argvals <- as.numeric(karcher.result$gammas$argvals)

  result <- .Call("wrap__alignment_phase_boxplot_rust",
                  gammas, argvals, as.numeric(factor))

  if (is.null(result)) {
    stop("phase.boxplot failed: check karcher.result contents")
  }

  structure(
    list(
      median = result$median,
      median.index = result$median_index,
      central.lower = result$central_lower,
      central.upper = result$central_upper,
      whisker.lower = result$whisker_lower,
      whisker.upper = result$whisker_upper,
      depths = result$depths,
      outlier.indices = result$outlier_indices,
      factor = result$factor,
      argvals = argvals,
      call = match.call()
    ),
    class = "phase.boxplot"
  )
}

#' Print Phase Boxplot
#'
#' @param x An object of class \code{phase.boxplot}.
#' @param ... Additional arguments (ignored).
#' @export
print.phase.boxplot <- function(x, ...) {
  cat("Phase Boxplot (Warping Functions)\n")
  cat("  Grid points:", length(x$median), "\n")
  cat("  Curves:", length(x$depths), "\n")
  cat("  Median index:", x$median.index, "\n")
  cat("  Outliers:", length(x$outlier.indices), "\n")
  cat("  Factor:", x$factor, "\n")
  invisible(x)
}

#' Plot Phase Boxplot
#'
#' Displays a functional boxplot of warping functions: median curve,
#' central 50% region (shaded), whiskers, and outliers in red.
#'
#' @param x An object of class \code{phase.boxplot}.
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object.
#'
#' @export
plot.phase.boxplot <- function(x, ...) {
  argvals <- x$argvals
  m <- length(argvals)

  df_band <- data.frame(
    t = argvals,
    median = x$median,
    central_lower = x$central.lower,
    central_upper = x$central.upper,
    whisker_lower = x$whisker.lower,
    whisker_upper = x$whisker.upper
  )

  p <- ggplot2::ggplot(df_band, ggplot2::aes(x = .data$t)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$whisker_lower, ymax = .data$whisker_upper),
      fill = "#4A90D9", alpha = 0.15
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$central_lower, ymax = .data$central_upper),
      fill = "#4A90D9", alpha = 0.35
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$median),
      color = "black", linewidth = 1.2
    )

  # Add outlier warps if present
  if (length(x$outlier.indices) > 0) {
    # Retrieve gammas from parent karcher.result is not feasible,
    # so we plot outlier indices as annotation
    # Instead, note: outlier.indices reference the original gamma matrix
  }

  p <- p +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "grey50") +
    ggplot2::labs(
      title = "Phase Boxplot",
      subtitle = paste0("Central 50% (dark), whiskers (light), ",
                        length(x$outlier.indices), " outlier(s)"),
      x = "t", y = expression(gamma(t))
    ) +
    ggplot2::theme_minimal()

  p
}

# =============================================================================
# Phase 7: Warp Inversion, Penalized Alignment, Diagnostics
# =============================================================================

#' Invert a Warping Function
#'
#' Compute the functional inverse of a warping function \eqn{\gamma},
#' i.e., find \eqn{\gamma^{-1}} such that
#' \eqn{\gamma(\gamma^{-1}(t)) = t}.
#'
#' @param gamma Numeric vector of warping function values.
#' @param argvals Numeric vector of evaluation points.
#'
#' @return Numeric vector of the inverse warping function.
#'
#' @seealso \code{\link{warp.inverse.error}}, \code{\link{warp.compose}}
#'
#' @export
#' @examples
#' t <- seq(0, 1, length.out = 50)
#' gamma <- t^2
#' gamma_inv <- warp.inverse(gamma, t)
warp.inverse <- function(gamma, argvals) {
  gamma <- as.numeric(gamma)
  argvals <- as.numeric(argvals)

  result <- .Call("wrap__alignment_invert_warp_rust", gamma, argvals)

  if (is.null(result)) {
    stop("warp.inverse failed: check gamma and argvals dimensions")
  }

  as.numeric(result)
}

#' Warp Inverse Error
#'
#' Compute the maximum reconstruction error of a warp inversion:
#' \eqn{\max_t |\gamma(\gamma^{-1}(t)) - t|}.
#'
#' @param gamma Numeric vector of warping function values.
#' @param argvals Numeric vector of evaluation points.
#'
#' @return A scalar error value.
#'
#' @seealso \code{\link{warp.inverse}}
#'
#' @export
#' @examples
#' t <- seq(0, 1, length.out = 50)
#' gamma <- t^2
#' warp.inverse.error(gamma, t)
warp.inverse.error <- function(gamma, argvals) {
  gamma <- as.numeric(gamma)
  argvals <- as.numeric(argvals)

  result <- .Call("wrap__alignment_warp_inverse_error_rust", gamma, argvals)
  as.numeric(result)
}

#' Penalized Elastic Pairwise Alignment
#'
#' Aligns curve \code{f2} to curve \code{f1} using elastic alignment with
#' a configurable roughness penalty on the warping function.
#'
#' @param f1 Numeric vector (reference curve).
#' @param f2 Numeric vector (curve to align).
#' @param argvals Numeric vector of evaluation points. If \code{NULL},
#'   defaults to \code{seq(0, 1, length.out = length(f1))}.
#' @param penalty Character: penalty type. One of \code{"first_order"},
#'   \code{"second_order"}, or \code{"combined"}.
#' @param lambda Numeric penalty weight (default 1.0).
#'
#' @return A list with components:
#' \describe{
#'   \item{gamma}{Warping function}
#'   \item{f.aligned}{Aligned version of f2}
#'   \item{distance}{Elastic distance after alignment}
#' }
#'
#' @seealso \code{\link{elastic.pair}}, \code{\link{elastic.lambda.cv}}
#'
#' @export
#' @examples
#' \donttest{
#' t <- seq(0, 1, length.out = 50)
#' f1 <- sin(2 * pi * t)
#' f2 <- sin(2 * pi * t^1.5)
#' res <- elastic.pair.penalized(f1, f2, argvals = t)
#' }
elastic.pair.penalized <- function(f1, f2, argvals = NULL,
                                    penalty = c("first_order", "second_order",
                                                "combined"),
                                    lambda = 1.0) {
  f1 <- as.numeric(f1)
  f2 <- as.numeric(f2)
  penalty <- match.arg(penalty)

  if (is.null(argvals)) {
    argvals <- seq(0, 1, length.out = length(f1))
  }
  argvals <- as.numeric(argvals)

  result <- .Call("wrap__alignment_elastic_pair_penalized_rust",
                  f1, f2, argvals, penalty, as.numeric(lambda))

  list(
    gamma = result$gamma,
    f.aligned = result$f_aligned,
    distance = result$distance
  )
}

#' Alignment Diagnostics
#'
#' Compute per-curve diagnostic information from a Karcher mean alignment,
#' identifying curves that may be poorly aligned (under-aligned, over-aligned,
#' or with non-monotone warps).
#'
#' @param fdataobj An object of class \code{fdata} (original data).
#' @param karcher.result An object of class \code{karcher.mean}.
#'
#' @return An object of class \code{alignment.diagnostics} with components:
#' \describe{
#'   \item{diagnostics}{List of per-curve diagnostic records, each containing
#'     \code{curve.index}, \code{warp.complexity}, \code{warp.smoothness},
#'     \code{is.under.aligned}, \code{is.over.aligned},
#'     \code{has.non.monotone}, \code{residual}, \code{distance.ratio},
#'     \code{flagged}, and \code{issues}}
#'   \item{flagged.indices}{Integer vector of 1-indexed flagged curve indices}
#'   \item{n.flagged}{Number of flagged curves}
#'   \item{health.score}{Overall alignment health score in \eqn{[0, 1]}}
#'   \item{call}{The matched call}
#' }
#'
#' @seealso \code{\link{karcher.mean}}, \code{\link{alignment.quality}},
#'   \code{\link{alignment.diagnostics.pairwise}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 15, 50)
#' for (i in 1:15) X[i, ] <- sin(2 * pi * (t - i / 60))
#' fd <- fdata(X, argvals = t)
#' km <- karcher.mean(fd, max.iter = 5)
#' diag <- alignment.diagnostics(fd, km)
#' diag
#' }
alignment.diagnostics <- function(fdataobj, karcher.result) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }
  if (!inherits(karcher.result, "karcher.mean")) {
    stop("karcher.result must be of class 'karcher.mean'")
  }

  argvals <- as.numeric(fdataobj$argvals)
  km_mean <- as.numeric(karcher.result$mean$data[1, ])
  km_mean_srsf <- karcher.result$mean_srsf
  if (is.null(km_mean_srsf)) {
    km_mean_srsf <- as.numeric(
      alignment_srsf_transform(karcher.result$mean$data, argvals)[1, ]
    )
  }

  result <- .Call("wrap__alignment_diagnose_rust",
                  fdataobj$data,
                  km_mean,
                  as.numeric(km_mean_srsf),
                  karcher.result$gammas$data,
                  karcher.result$aligned$data,
                  as.integer(karcher.result$n.iter),
                  as.logical(karcher.result$converged),
                  argvals)

  if (is.null(result)) {
    stop("alignment.diagnostics failed: check data dimensions")
  }

  structure(
    list(
      diagnostics = result$diagnostics,
      flagged.indices = result$flagged_indices,
      n.flagged = result$n_flagged,
      health.score = result$health_score,
      call = match.call()
    ),
    class = "alignment.diagnostics"
  )
}

#' Print Alignment Diagnostics
#'
#' @param x An object of class \code{alignment.diagnostics}.
#' @param ... Additional arguments (ignored).
#' @export
print.alignment.diagnostics <- function(x, ...) {
  n <- length(x$diagnostics)
  cat("Alignment Diagnostics\n")
  cat("  Curves:", n, "\n")
  cat("  Flagged:", x$n.flagged, "of", n,
      paste0("(", format(100 * x$n.flagged / n, digits = 3), "%)"), "\n")
  cat("  Health score:", format(x$health.score, digits = 4), "\n")
  if (x$n.flagged > 0) {
    cat("  Flagged indices:", paste(x$flagged.indices, collapse = ", "), "\n")
  }
  invisible(x)
}

#' Plot Alignment Diagnostics
#'
#' Displays a per-curve diagnostic summary showing the distance ratio
#' (before/after alignment variance ratio) for each curve, with flagged
#' curves highlighted.
#'
#' @param x An object of class \code{alignment.diagnostics}.
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object.
#'
#' @export
plot.alignment.diagnostics <- function(x, ...) {
  n <- length(x$diagnostics)
  curve_idx <- integer(n)
  dist_ratio <- numeric(n)
  flagged <- logical(n)

  for (i in seq_len(n)) {
    d <- x$diagnostics[[i]]
    curve_idx[i] <- d$curve_index
    dist_ratio[i] <- d$distance_ratio
    flagged[i] <- as.logical(d$flagged)
  }

  df <- data.frame(
    curve = curve_idx,
    distance_ratio = dist_ratio,
    flagged = flagged
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$curve,
                                         y = .data$distance_ratio)) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = .data$curve, y = 0, yend = .data$distance_ratio,
                   color = .data$flagged),
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$flagged), size = 2
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "#4A90D9", "TRUE" = "#D55E00"),
      labels = c("FALSE" = "OK", "TRUE" = "Flagged"),
      name = "Status"
    ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    ggplot2::labs(
      title = paste0("Alignment Diagnostics (health = ",
                     format(x$health.score, digits = 3), ")"),
      x = "Curve Index", y = "Distance Ratio"
    ) +
    ggplot2::theme_minimal()

  p
}

#' Pairwise Alignment Diagnostics
#'
#' Compute diagnostic information for a single pairwise alignment result.
#'
#' @param f1 Numeric vector (reference curve).
#' @param f2 Numeric vector (curve that was aligned).
#' @param result A pairwise alignment result list with components
#'   \code{$gamma} (or \code{$gamma}), \code{$f.aligned} (or
#'   \code{$f_aligned}), and \code{$distance}.
#' @param argvals Numeric vector of evaluation points.
#'
#' @return A list with diagnostic fields:
#'   \code{curve.index}, \code{warp.complexity}, \code{warp.smoothness},
#'   \code{is.under.aligned}, \code{is.over.aligned},
#'   \code{has.non.monotone}, \code{residual}, \code{distance.ratio},
#'   \code{flagged}, \code{issues}.
#'
#' @seealso \code{\link{elastic.pair}}, \code{\link{alignment.diagnostics}}
#'
#' @export
#' @examples
#' \donttest{
#' t <- seq(0, 1, length.out = 50)
#' f1 <- sin(2 * pi * t)
#' f2 <- sin(2 * pi * t^1.5)
#' fd1 <- fdata(matrix(f1, 1, 50), argvals = t)
#' fd2 <- fdata(matrix(f2, 1, 50), argvals = t)
#' res <- elastic.pair(fd1, fd2)
#' diag <- alignment.diagnostics.pairwise(f1, f2, res, t)
#' }
alignment.diagnostics.pairwise <- function(f1, f2, result, argvals) {
  f1 <- as.numeric(f1)
  f2 <- as.numeric(f2)
  argvals <- as.numeric(argvals)

  # Accept both dot and underscore naming from result
  gamma <- if (!is.null(result$gamma)) result$gamma else result$gamma
  f_aligned <- if (!is.null(result$f.aligned)) {
    as.numeric(result$f.aligned$data[1, ])
  } else if (!is.null(result$f_aligned)) {
    as.numeric(result$f_aligned)
  } else {
    stop("result must have $f.aligned or $f_aligned")
  }
  distance <- result$distance

  res <- .Call("wrap__alignment_diagnose_pairwise_rust",
               f1, f2,
               as.numeric(gamma), as.numeric(f_aligned),
               as.numeric(distance), argvals)

  if (is.null(res)) {
    stop("alignment.diagnostics.pairwise failed: check dimensions")
  }

  list(
    curve.index = res$curve_index,
    warp.complexity = res$warp_complexity,
    warp.smoothness = res$warp_smoothness,
    is.under.aligned = as.logical(res$is_under_aligned),
    is.over.aligned = as.logical(res$is_over_aligned),
    has.non.monotone = as.logical(res$has_non_monotone),
    residual = res$residual,
    distance.ratio = res$distance_ratio,
    flagged = as.logical(res$flagged),
    issues = res$issues
  )
}

# =============================================================================
# Phase 8: Shape Analysis
# =============================================================================

#' Shape Representative (Orbit Representative)
#'
#' Compute the canonical representative of a curve in a quotient space.
#' Depending on the quotient, this factors out reparameterization,
#' translation, and/or scale.
#'
#' @param f Numeric vector of curve values.
#' @param argvals Numeric vector of evaluation points. If \code{NULL},
#'   defaults to \code{seq(0, 1, length.out = length(f))}.
#' @param quotient Character: the quotient space. One of
#'   \code{"reparameterization"}, \code{"translation"}, or \code{"scale"}.
#'
#' @return A list with components:
#' \describe{
#'   \item{representative}{The canonical representative curve}
#'   \item{representative.srsf}{SRSF of the representative}
#'   \item{gamma}{Warping function to the representative}
#'   \item{translation}{Translation applied (0 for reparameterization-only)}
#'   \item{scale}{Scale factor applied (1 for reparameterization-only)}
#' }
#'
#' @references
#' Srivastava, A. and Klassen, E. (2016). \emph{Functional and Shape Data
#' Analysis}. Springer.
#'
#' @seealso \code{\link{shape.distance}}, \code{\link{shape.mean}}
#'
#' @export
#' @examples
#' \donttest{
#' t <- seq(0, 1, length.out = 50)
#' f <- sin(2 * pi * t)
#' rep <- shape.representative(f, argvals = t)
#' }
shape.representative <- function(f, argvals = NULL,
                                  quotient = c("reparameterization",
                                               "translation", "scale")) {
  f <- as.numeric(f)
  quotient <- match.arg(quotient)

  if (is.null(argvals)) {
    argvals <- seq(0, 1, length.out = length(f))
  }
  argvals <- as.numeric(argvals)

  result <- .Call("wrap__alignment_orbit_representative_rust",
                  f, argvals, quotient)

  if (is.null(result)) {
    stop("shape.representative failed: check input dimensions")
  }

  list(
    representative = result$representative,
    representative.srsf = result$representative_srsf,
    gamma = result$gamma,
    translation = result$translation,
    scale = result$scale
  )
}

#' Shape Distance Between Two Curves
#'
#' Compute the elastic shape distance between two curves in a specified
#' quotient space. The shape distance factors out the specified
#' nuisance transformations (reparameterization, translation, scale).
#'
#' @param f1 Numeric vector (first curve).
#' @param f2 Numeric vector (second curve).
#' @param argvals Numeric vector of evaluation points. If \code{NULL},
#'   defaults to \code{seq(0, 1, length.out = length(f1))}.
#' @param quotient Character: the quotient space. One of
#'   \code{"reparameterization"}, \code{"translation"}, or \code{"scale"}.
#' @param lambda Regularization parameter (default 0).
#'
#' @return A list with components:
#' \describe{
#'   \item{distance}{Shape distance}
#'   \item{gamma}{Optimal warping function}
#'   \item{f2.aligned}{Aligned version of f2}
#' }
#'
#' @references
#' Srivastava, A. and Klassen, E. (2016). \emph{Functional and Shape Data
#' Analysis}. Springer.
#'
#' @seealso \code{\link{shape.representative}}, \code{\link{shape.mean}},
#'   \code{\link{shape.distance.matrix}}
#'
#' @export
#' @examples
#' \donttest{
#' t <- seq(0, 1, length.out = 50)
#' f1 <- sin(2 * pi * t)
#' f2 <- 2 * sin(2 * pi * t^1.5) + 1
#' d <- shape.distance(f1, f2, argvals = t, quotient = "scale")
#' d$distance
#' }
shape.distance <- function(f1, f2, argvals = NULL,
                            quotient = c("reparameterization",
                                         "translation", "scale"),
                            lambda = 0) {
  f1 <- as.numeric(f1)
  f2 <- as.numeric(f2)
  quotient <- match.arg(quotient)

  if (is.null(argvals)) {
    argvals <- seq(0, 1, length.out = length(f1))
  }
  argvals <- as.numeric(argvals)

  result <- .Call("wrap__alignment_shape_distance_rust",
                  f1, f2, argvals, quotient, as.numeric(lambda))

  if (is.null(result)) {
    stop("shape.distance failed: check input dimensions")
  }

  list(
    distance = result$distance,
    gamma = result$gamma,
    f2.aligned = result$f2_aligned
  )
}

#' Shape Mean (Karcher Mean in Quotient Space)
#'
#' Compute the Karcher (Frechet) mean of functional data in a shape
#' quotient space. This simultaneously estimates the mean shape and
#' aligns all curves, factoring out the specified nuisance transformations.
#'
#' @param fdataobj An object of class \code{fdata}.
#' @param quotient Character: the quotient space. One of
#'   \code{"reparameterization"}, \code{"translation"}, or \code{"scale"}.
#' @param lambda Regularization parameter (default 0).
#' @param max.iter Maximum number of iterations (default 20).
#' @param tol Convergence tolerance (default 1e-4).
#'
#' @return An object of class \code{shape.mean} with components:
#' \describe{
#'   \item{mean}{The shape mean curve (numeric vector)}
#'   \item{mean.srsf}{SRSF of the mean curve}
#'   \item{gammas}{Matrix of warping functions (n x m)}
#'   \item{aligned.data}{Matrix of aligned curves (n x m)}
#'   \item{n.iter}{Number of iterations}
#'   \item{converged}{Logical: did the algorithm converge?}
#'   \item{fdataobj}{Original fdata object}
#'   \item{quotient}{Quotient space used}
#'   \item{call}{The matched call}
#' }
#'
#' @references
#' Srivastava, A. and Klassen, E. (2016). \emph{Functional and Shape Data
#' Analysis}. Springer.
#'
#' @seealso \code{\link{shape.distance}}, \code{\link{shape.representative}},
#'   \code{\link{shape.distance.matrix}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 10, 50)
#' for (i in 1:10) X[i, ] <- sin(2 * pi * (t - i / 50)) + rnorm(50, 0, 0.1)
#' fd <- fdata(X, argvals = t)
#' sm <- shape.mean(fd, quotient = "reparameterization", max.iter = 10)
#' sm
#' }
shape.mean <- function(fdataobj,
                        quotient = c("reparameterization", "translation",
                                     "scale"),
                        lambda = 0, max.iter = 20, tol = 1e-4) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }
  quotient <- match.arg(quotient)

  argvals <- as.numeric(fdataobj$argvals)

  result <- .Call("wrap__alignment_shape_mean_rust",
                  fdataobj$data, argvals, quotient,
                  as.numeric(lambda), as.integer(max.iter), as.numeric(tol))

  if (is.null(result)) {
    stop("shape.mean failed: check data dimensions and parameters")
  }

  structure(
    list(
      mean = result$mean,
      mean.srsf = result$mean_srsf,
      gammas = result$gammas,
      aligned.data = result$aligned_data,
      n.iter = result$n_iter,
      converged = result$converged,
      fdataobj = fdataobj,
      quotient = quotient,
      call = match.call()
    ),
    class = "shape.mean"
  )
}

#' Print Shape Mean
#'
#' @param x An object of class \code{shape.mean}.
#' @param ... Additional arguments (ignored).
#' @export
print.shape.mean <- function(x, ...) {
  n <- nrow(x$fdataobj$data)
  m <- ncol(x$fdataobj$data)
  cat("Shape Mean (Quotient Space)\n")
  cat("  Curves:", n, "x", m, "grid points\n")
  cat("  Quotient:", x$quotient, "\n")
  cat("  Iterations:", x$n.iter, "\n")
  cat("  Converged:", x$converged, "\n")
  invisible(x)
}

#' Plot Shape Mean
#'
#' Displays the shape mean curve with aligned input curves in the background.
#'
#' @param x An object of class \code{shape.mean}.
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object.
#'
#' @export
plot.shape.mean <- function(x, ...) {
  argvals <- as.numeric(x$fdataobj$argvals)
  m <- length(argvals)

  # aligned.data is a matrix (n x m)
  aligned <- x$aligned.data
  if (is.matrix(aligned)) {
    n <- nrow(aligned)
  } else {
    n <- 1
    aligned <- matrix(aligned, nrow = 1)
  }

  df_aligned <- data.frame(
    curve_id = rep(seq_len(n), each = m),
    argval = rep(argvals, n),
    value = as.vector(t(aligned))
  )

  df_mean <- data.frame(
    argval = argvals,
    value = as.numeric(x$mean)
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
    ggplot2::labs(
      title = paste0("Shape Mean (", x$quotient, ")"),
      x = "t", y = "value"
    ) +
    ggplot2::theme_minimal()

  p
}

#' Shape Distance Matrix
#'
#' Compute the pairwise shape distance matrix for functional data in a
#' specified quotient space.
#'
#' @param fdataobj An object of class \code{fdata}.
#' @param quotient Character: the quotient space. One of
#'   \code{"reparameterization"}, \code{"translation"}, or \code{"scale"}.
#' @param lambda Regularization parameter (default 0).
#'
#' @return A symmetric distance matrix (n x n).
#'
#' @references
#' Srivastava, A. and Klassen, E. (2016). \emph{Functional and Shape Data
#' Analysis}. Springer.
#'
#' @seealso \code{\link{shape.distance}}, \code{\link{shape.mean}},
#'   \code{\link{elastic.distance}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 30)
#' X <- matrix(0, 8, 30)
#' for (i in 1:8) X[i, ] <- sin(2 * pi * (t - i / 40))
#' fd <- fdata(X, argvals = t)
#' D <- shape.distance.matrix(fd, quotient = "reparameterization")
#' }
shape.distance.matrix <- function(fdataobj,
                                   quotient = c("reparameterization",
                                                "translation", "scale"),
                                   lambda = 0) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }
  quotient <- match.arg(quotient)

  argvals <- as.numeric(fdataobj$argvals)

  D <- .Call("wrap__alignment_shape_self_distance_matrix_rust",
             fdataobj$data, argvals, quotient, as.numeric(lambda))

  D <- as.matrix(D)

  if (!is.null(rownames(fdataobj$data))) {
    rownames(D) <- rownames(fdataobj$data)
    colnames(D) <- rownames(fdataobj$data)
  }

  D
}

# =============================================================================
# Karcher Median
# =============================================================================

#' Karcher Median in Elastic Metric
#'
#' Compute the Karcher median of functional data in the elastic metric using the
#' Weiszfeld algorithm. The Karcher median minimises the sum of distances rather
#' than the sum of squared distances, making it more robust to outliers than the
#' Karcher mean.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param max.iter Maximum number of Weiszfeld iterations (default 30).
#' @param tol Convergence tolerance (default 1e-4).
#' @param lambda Regularisation parameter controlling warping smoothness
#'   (default 0).
#' @param trim Fraction of most distant curves to trim before aggregation
#'   (default 0, no trimming).
#'
#' @return An object of class 'karcher.median' with components:
#' \describe{
#'   \item{median}{fdata of the Karcher median curve}
#'   \item{median_srsf}{numeric vector of the median SRSF}
#'   \item{aligned}{fdata of aligned curves}
#'   \item{gammas}{fdata of warping functions}
#'   \item{n.iter}{number of iterations used}
#'   \item{converged}{logical indicating convergence}
#'   \item{fdataobj}{the original fdata input}
#' }
#'
#' @references
#' Fletcher, P.T., Venkatasubramanian, S., and Joshi, S. (2009). The geometric
#' median on Riemannian manifolds with application to robust atlas estimation.
#' \emph{NeuroImage}, 45(1):S143--S152.
#'
#' Srivastava, A. and Klassen, E. (2016). \emph{Functional and Shape Data
#' Analysis}. Springer.
#'
#' @seealso \code{\link{karcher.mean}}, \code{\link{robust.karcher.mean}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 15, 50)
#' for (i in 1:15) X[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
#' fd <- fdata(X, argvals = t)
#' km <- karcher.median(fd)
#' }
karcher.median <- function(fdataobj, max.iter = 30, tol = 1e-4,
                           lambda = 0, trim = 0) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  argvals <- fdataobj$argvals
  res <- alignment_karcher_median(fdataobj$data, argvals,
                                  as.integer(max.iter), as.numeric(tol),
                                  as.numeric(lambda), as.numeric(trim))

  result <- list(
    median = fdata(matrix(res$mean, nrow = 1), argvals = argvals),
    median_srsf = res$mean_srsf,
    aligned = fdata(res$aligned_data, argvals = argvals),
    gammas = fdata(res$gammas, argvals = argvals),
    n.iter = res$n_iter,
    converged = res$converged,
    fdataobj = fdataobj
  )
  class(result) <- "karcher.median"
  result
}

#' @export
print.karcher.median <- function(x, ...) {
  n <- nrow(x$fdataobj$data)
  m <- ncol(x$fdataobj$data)
  cat("Karcher Median (Elastic)\n")
  cat("  Curves:", n, "x", m, "grid points\n")
  cat("  Iterations:", x$n.iter, "\n")
  cat("  Converged:", x$converged, "\n")
  invisible(x)
}

# =============================================================================
# Robust Karcher Mean
# =============================================================================

#' Trimmed Karcher Mean in Elastic Metric
#'
#' Compute a robust version of the Karcher mean by iteratively trimming the most
#' distant curves before updating the mean estimate. This reduces the influence
#' of outlying observations on the mean shape.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param trim Fraction of most distant curves to trim at each iteration
#'   (default 0.1).
#' @param max.iter Maximum number of iterations (default 20).
#' @param tol Convergence tolerance (default 1e-4).
#' @param lambda Regularisation parameter controlling warping smoothness
#'   (default 0).
#'
#' @return An object of class 'karcher.mean' with components:
#' \describe{
#'   \item{mean}{fdata of the trimmed Karcher mean curve}
#'   \item{mean_srsf}{numeric vector of the mean SRSF}
#'   \item{aligned}{fdata of aligned curves}
#'   \item{gammas}{fdata of warping functions}
#'   \item{n.iter}{number of iterations used}
#'   \item{converged}{logical indicating convergence}
#'   \item{fdataobj}{the original fdata input}
#' }
#'
#' @references
#' Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
#' functional data using phase and amplitude separation.
#' \emph{Computational Statistics & Data Analysis}, 61:50--66.
#'
#' @seealso \code{\link{karcher.mean}}, \code{\link{karcher.median}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 20, 50)
#' for (i in 1:20) X[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
#' fd <- fdata(X, argvals = t)
#' rkm <- robust.karcher.mean(fd, trim = 0.1)
#' }
robust.karcher.mean <- function(fdataobj, trim = 0.1, max.iter = 20,
                                tol = 1e-4, lambda = 0) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  argvals <- fdataobj$argvals
  res <- alignment_robust_karcher_mean(fdataobj$data, argvals,
                                       as.integer(max.iter), as.numeric(tol),
                                       as.numeric(lambda), as.numeric(trim))

  result <- list(
    mean = fdata(matrix(res$mean, nrow = 1), argvals = argvals),
    mean_srsf = res$mean_srsf,
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
# Elastic Depth
# =============================================================================

#' Elastic Depth for Functional Data
#'
#' Compute elastic depth scores that decompose into amplitude and phase
#' components. Elastic depth measures the centrality of each curve relative to
#' the sample, accounting for both vertical variation (amplitude) and horizontal
#' variation (phase/warping).
#'
#' @param fdataobj An object of class 'fdata'.
#' @param lambda Regularisation parameter controlling warping smoothness
#'   (default 0).
#'
#' @return An object of class 'elastic.depth' with components:
#' \describe{
#'   \item{amplitude_depth}{numeric vector of amplitude depth scores}
#'   \item{phase_depth}{numeric vector of phase depth scores}
#'   \item{combined_depth}{numeric vector of combined depth scores}
#' }
#'
#' @references
#' Lopez-Pintado, S. and Romo, J. (2009). On the concept of depth for
#' functional data. \emph{Journal of the American Statistical Association},
#' 104(486):718--734.
#'
#' Srivastava, A. and Klassen, E. (2016). \emph{Functional and Shape Data
#' Analysis}. Springer.
#'
#' @seealso \code{\link{elastic.outlier.detection}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 15, 50)
#' for (i in 1:15) X[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
#' fd <- fdata(X, argvals = t)
#' ed <- elastic.depth(fd)
#' }
elastic.depth <- function(fdataobj, lambda = 0) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  res <- alignment_elastic_depth(fdataobj$data, fdataobj$argvals,
                                 as.numeric(lambda))

  result <- list(
    amplitude_depth = res$amplitude_depth,
    phase_depth = res$phase_depth,
    combined_depth = res$combined_depth
  )
  class(result) <- "elastic.depth"
  result
}

#' @export
print.elastic.depth <- function(x, ...) {
  n <- length(x$combined_depth)
  cat("Elastic Depth\n")
  cat("  Curves:", n, "\n")
  cat("  Amplitude depth range:",
      round(min(x$amplitude_depth), 4), "-",
      round(max(x$amplitude_depth), 4), "\n")
  cat("  Phase depth range:",
      round(min(x$phase_depth), 4), "-",
      round(max(x$phase_depth), 4), "\n")
  invisible(x)
}

# =============================================================================
# Elastic Outlier Detection
# =============================================================================

#' Elastic Outlier Detection
#'
#' Identify outlying curves using SRVF-based elastic distances and Tukey's fence
#' rule. Outliers are detected in both the amplitude and phase domains.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param lambda Regularisation parameter controlling warping smoothness
#'   (default 0).
#' @param alpha Significance level for the Tukey fence (default 0.05).
#' @param use.median Logical; if TRUE, use the Karcher median instead of the
#'   Karcher mean as the central tendency reference (default FALSE).
#'
#' @return A list with components:
#' \describe{
#'   \item{outlier_indices}{integer vector of detected outlier indices}
#'   \item{amplitude_outliers}{integer vector of amplitude outlier indices}
#'   \item{phase_outliers}{integer vector of phase outlier indices}
#'   \item{amplitude_distances}{numeric vector of amplitude distances}
#'   \item{phase_distances}{numeric vector of phase distances}
#'   \item{amplitude_fence}{numeric; the amplitude fence threshold}
#'   \item{phase_fence}{numeric; the phase fence threshold}
#' }
#'
#' @references
#' Xie, W., Kurtek, S., Bharath, K., and Sun, Y. (2017). A geometric approach
#' to visualization of variability in functional data.
#' \emph{Journal of the American Statistical Association},
#' 112(519):979--993.
#'
#' @seealso \code{\link{elastic.depth}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 20, 50)
#' for (i in 1:20) X[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
#' # Add an outlier
#' X[20, ] <- 5 * cos(4 * pi * t)
#' fd <- fdata(X, argvals = t)
#' out <- elastic.outlier.detection(fd)
#' }
elastic.outlier.detection <- function(fdataobj, lambda = 0, alpha = 0.05,
                                      use.median = FALSE) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  res <- alignment_elastic_outlier(fdataobj$data, fdataobj$argvals,
                                   as.numeric(lambda), as.numeric(alpha),
                                   as.logical(use.median))
  res
}

# =============================================================================
# Shape Confidence Interval
# =============================================================================

#' Bootstrap Confidence Bands for Elastic Karcher Mean
#'
#' Construct pointwise bootstrap confidence bands for the elastic Karcher mean
#' by resampling curves and recomputing the mean for each bootstrap replicate.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param n.bootstrap Number of bootstrap replicates (default 200).
#' @param confidence.level Confidence level for the bands (default 0.95).
#' @param lambda Regularisation parameter controlling warping smoothness
#'   (default 0).
#' @param max.iter Maximum number of Karcher mean iterations per replicate
#'   (default 15).
#' @param tol Convergence tolerance for each Karcher mean computation
#'   (default 1e-4).
#' @param seed Random seed for reproducibility (default 42).
#'
#' @return An object of class 'shape.ci' with components:
#' \describe{
#'   \item{mean}{fdata of the Karcher mean curve}
#'   \item{lower}{numeric vector of lower confidence band values}
#'   \item{upper}{numeric vector of upper confidence band values}
#'   \item{confidence.level}{the confidence level used}
#'   \item{n.bootstrap}{the number of bootstrap replicates used}
#'   \item{argvals}{the evaluation grid}
#' }
#'
#' @references
#' Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
#' functional data using phase and amplitude separation.
#' \emph{Computational Statistics & Data Analysis}, 61:50--66.
#'
#' @seealso \code{\link{karcher.mean}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 20, 50)
#' for (i in 1:20) X[i, ] <- sin(2 * pi * t + runif(1, -0.2, 0.2))
#' fd <- fdata(X, argvals = t)
#' ci <- shape.confidence.interval(fd, n.bootstrap = 50)
#' }
shape.confidence.interval <- function(fdataobj, n.bootstrap = 200,
                                      confidence.level = 0.95, lambda = 0,
                                      max.iter = 15, tol = 1e-4, seed = 42) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  argvals <- fdataobj$argvals
  res <- alignment_shape_ci(fdataobj$data, argvals,
                            as.integer(n.bootstrap),
                            as.numeric(confidence.level),
                            as.numeric(lambda),
                            as.integer(max.iter),
                            as.numeric(tol),
                            as.integer(seed))

  result <- list(
    mean = fdata(matrix(res$mean, nrow = 1), argvals = argvals),
    lower = res$lower,
    upper = res$upper,
    confidence.level = confidence.level,
    n.bootstrap = n.bootstrap,
    argvals = argvals
  )
  class(result) <- "shape.ci"
  result
}

#' @export
print.shape.ci <- function(x, ...) {
  cat("Shape Confidence Interval\n")
  cat("  Confidence level:", x$confidence.level, "\n")
  cat("  Bootstrap replicates:", x$n.bootstrap, "\n")
  cat("  Grid points:", length(x$argvals), "\n")
  invisible(x)
}

# =============================================================================
# Bayesian Pairwise Alignment
# =============================================================================

#' Bayesian Pairwise Curve Alignment
#'
#' Align two curves using a Bayesian framework with MCMC sampling of the warping
#' function posterior. This provides uncertainty quantification for the alignment
#' via posterior samples of the warping function.
#'
#' @param f1 Numeric vector of the first curve's values.
#' @param f2 Numeric vector of the second curve's values.
#' @param argvals Numeric vector of the evaluation grid (common to both curves).
#' @param n.samples Number of MCMC samples to draw (default 2000).
#' @param burn.in Number of initial samples to discard as burn-in (default 500).
#' @param step.size Step size for the MCMC proposal (default 0.1).
#' @param proposal.variance Variance of the MCMC proposal distribution
#'   (default 1.0).
#' @param seed Random seed for reproducibility (default 42).
#'
#' @return A list with components:
#' \describe{
#'   \item{gamma.mean}{numeric vector of the posterior mean warping function}
#'   \item{gamma.samples}{matrix of posterior warping function samples
#'     (samples x grid)}
#'   \item{f2.aligned}{numeric vector of the aligned second curve}
#'   \item{acceptance.rate}{the MCMC acceptance rate}
#' }
#'
#' @references
#' Cheng, W., Dryden, I.L., and Huang, X. (2016). Bayesian registration of
#' functions and curves. \emph{Bayesian Analysis}, 11(2):447--475.
#'
#' @seealso \code{\link{elastic.align}}
#'
#' @export
#' @examples
#' \donttest{
#' t <- seq(0, 1, length.out = 100)
#' f1 <- sin(2 * pi * t)
#' f2 <- sin(2 * pi * (t^1.3))
#' res <- bayesian.align.pair(f1, f2, t, n.samples = 500, burn.in = 100)
#' }
bayesian.align.pair <- function(f1, f2, argvals, n.samples = 2000,
                                burn.in = 500, step.size = 0.1,
                                proposal.variance = 1.0, seed = 42) {
  res <- alignment_bayesian_pair(as.numeric(f1), as.numeric(f2),
                                 as.numeric(argvals),
                                 as.integer(n.samples),
                                 as.integer(burn.in),
                                 as.numeric(step.size),
                                 as.numeric(proposal.variance),
                                 as.integer(seed))
  res
}

# =============================================================================
# Multi-resolution Pairwise Alignment
# =============================================================================

#' Multi-Resolution Pairwise Curve Alignment
#'
#' Align two curves using a coarse-to-fine multi-resolution strategy. The curves
#' are first coarsened, aligned at low resolution, and then the alignment is
#' progressively refined to the original grid resolution. This can improve
#' robustness and speed for highly misaligned curves.
#'
#' @param f1 Numeric vector of the first curve's values.
#' @param f2 Numeric vector of the second curve's values.
#' @param argvals Numeric vector of the evaluation grid (common to both curves).
#' @param coarsen.factor Coarsening factor for the initial resolution
#'   (default 4).
#' @param n.refine Number of refinement steps (default 10).
#' @param step.size Step size for gradient-based refinement (default 0.01).
#' @param lambda Regularisation parameter controlling warping smoothness
#'   (default 0).
#'
#' @return A list with components:
#' \describe{
#'   \item{gamma}{numeric vector of the estimated warping function}
#'   \item{f2.aligned}{numeric vector of the aligned second curve}
#'   \item{distance}{elastic distance after alignment}
#' }
#'
#' @seealso \code{\link{elastic.align}}
#'
#' @export
#' @examples
#' \donttest{
#' t <- seq(0, 1, length.out = 200)
#' f1 <- sin(2 * pi * t)
#' f2 <- sin(2 * pi * (t^1.5))
#' res <- elastic.align.pair.multires(f1, f2, t)
#' }
elastic.align.pair.multires <- function(f1, f2, argvals, coarsen.factor = 4,
                                        n.refine = 10, step.size = 0.01,
                                        lambda = 0) {
  res <- alignment_multires_pair(as.numeric(f1), as.numeric(f2),
                                 as.numeric(argvals),
                                 as.integer(coarsen.factor),
                                 as.integer(n.refine),
                                 as.numeric(step.size),
                                 as.numeric(lambda))
  res
}

# =============================================================================
# Closed Curve Alignment
# =============================================================================

#' Elastic Alignment of Closed Curves
#'
#' Align two closed curves in the elastic framework. Closed curves satisfy
#' \eqn{f(0) = f(1)} and the alignment optimises over both reparameterisation
#' and rotation (circular shift) of the parameter domain.
#'
#' @param f1 Numeric vector of the first closed curve's values.
#' @param f2 Numeric vector of the second closed curve's values.
#' @param argvals Numeric vector of the evaluation grid (common to both curves).
#' @param lambda Regularisation parameter controlling warping smoothness
#'   (default 0).
#'
#' @return A list with components:
#' \describe{
#'   \item{gamma}{numeric vector of the estimated warping function}
#'   \item{f2.aligned}{numeric vector of the aligned second curve}
#'   \item{distance}{elastic distance after alignment}
#'   \item{rotation}{the optimal rotation applied}
#' }
#'
#' @references
#' Srivastava, A. and Klassen, E. (2016). \emph{Functional and Shape Data
#' Analysis}. Springer.
#'
#' @seealso \code{\link{elastic.distance.closed}},
#'   \code{\link{karcher.mean.closed}}
#'
#' @export
#' @examples
#' \donttest{
#' t <- seq(0, 1, length.out = 100)
#' f1 <- sin(2 * pi * t)
#' f2 <- sin(2 * pi * t + 0.5)
#' res <- elastic.align.pair.closed(f1, f2, t)
#' }
elastic.align.pair.closed <- function(f1, f2, argvals, lambda = 0) {
  res <- alignment_closed_pair(as.numeric(f1), as.numeric(f2),
                               as.numeric(argvals), as.numeric(lambda))
  res
}

# =============================================================================
# Closed Curve Distance
# =============================================================================

#' Elastic Distance Between Closed Curves
#'
#' Compute the elastic (Fisher-Rao) distance between two closed curves,
#' optimising over both reparameterisation and rotation.
#'
#' @param f1 Numeric vector of the first closed curve's values.
#' @param f2 Numeric vector of the second closed curve's values.
#' @param argvals Numeric vector of the evaluation grid (common to both curves).
#' @param lambda Regularisation parameter controlling warping smoothness
#'   (default 0).
#'
#' @return Numeric scalar; the elastic distance between the two closed curves.
#'
#' @references
#' Srivastava, A. and Klassen, E. (2016). \emph{Functional and Shape Data
#' Analysis}. Springer.
#'
#' @seealso \code{\link{elastic.align.pair.closed}},
#'   \code{\link{karcher.mean.closed}}
#'
#' @export
#' @examples
#' \donttest{
#' t <- seq(0, 1, length.out = 100)
#' f1 <- sin(2 * pi * t)
#' f2 <- sin(2 * pi * t + 0.5)
#' d <- elastic.distance.closed(f1, f2, t)
#' }
elastic.distance.closed <- function(f1, f2, argvals, lambda = 0) {
  alignment_closed_distance(as.numeric(f1), as.numeric(f2),
                            as.numeric(argvals), as.numeric(lambda))
}

# =============================================================================
# Karcher Mean for Closed Curves
# =============================================================================

#' Karcher Mean for Closed Curves
#'
#' Compute the Karcher (Frechet) mean for a sample of closed curves in the
#' elastic metric. This optimises simultaneously over reparameterisations and
#' rotations of the parameter domain.
#'
#' @param fdataobj An object of class 'fdata' containing closed curves.
#' @param max.iter Maximum number of iterations (default 15).
#' @param tol Convergence tolerance (default 1e-4).
#' @param lambda Regularisation parameter controlling warping smoothness
#'   (default 0).
#'
#' @return An object of class 'karcher.mean' with components:
#' \describe{
#'   \item{mean}{fdata of the Karcher mean curve}
#'   \item{mean_srsf}{numeric vector of the mean SRSF}
#'   \item{aligned}{fdata of aligned curves}
#'   \item{gammas}{fdata of warping functions}
#'   \item{n.iter}{number of iterations used}
#'   \item{converged}{logical indicating convergence}
#'   \item{fdataobj}{the original fdata input}
#' }
#'
#' @references
#' Srivastava, A. and Klassen, E. (2016). \emph{Functional and Shape Data
#' Analysis}. Springer.
#'
#' @seealso \code{\link{karcher.mean}}, \code{\link{elastic.align.pair.closed}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 100)
#' X <- matrix(0, 10, 100)
#' for (i in 1:10) X[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
#' fd <- fdata(X, argvals = t)
#' km <- karcher.mean.closed(fd)
#' }
karcher.mean.closed <- function(fdataobj, max.iter = 15, tol = 1e-4,
                                lambda = 0) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  argvals <- fdataobj$argvals
  res <- alignment_karcher_mean_closed(fdataobj$data, argvals,
                                       as.integer(max.iter), as.numeric(tol),
                                       as.numeric(lambda))

  result <- list(
    mean = fdata(matrix(res$mean, nrow = 1), argvals = argvals),
    mean_srsf = res$mean_srsf,
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
# Elastic Partial Matching
# =============================================================================

#' Elastic Partial Curve Matching
#'
#' Find the optimal partial match of a template curve within a longer target
#' curve using the elastic framework. This identifies the sub-interval of the
#' target that best matches the template, along with the optimal warping.
#'
#' @param template Numeric vector of the template curve's values.
#' @param target Numeric vector of the target curve's values.
#' @param argvals.template Numeric vector of evaluation points for the template.
#' @param argvals.target Numeric vector of evaluation points for the target.
#' @param lambda Regularisation parameter controlling warping smoothness
#'   (default 0).
#' @param min.span Minimum fraction of the target domain that the matched
#'   sub-interval must span (default 0.3).
#'
#' @return A list with components:
#' \describe{
#'   \item{gamma}{numeric vector of the warping function on the matched
#'     sub-interval}
#'   \item{match_start}{start of the matched sub-interval on the target domain}
#'   \item{match_end}{end of the matched sub-interval on the target domain}
#'   \item{distance}{elastic distance of the partial match}
#'   \item{target_aligned}{numeric vector of the matched and aligned target
#'     segment}
#' }
#'
#' @references
#' Srivastava, A. and Klassen, E. (2016). \emph{Functional and Shape Data
#' Analysis}. Springer.
#'
#' @export
#' @examples
#' \donttest{
#' t1 <- seq(0, 1, length.out = 50)
#' t2 <- seq(0, 2, length.out = 100)
#' template <- sin(2 * pi * t1)
#' target <- c(rep(0, 20), sin(2 * pi * seq(0, 1, length.out = 50)), rep(0, 30))
#' res <- elastic.partial.match(template, target, t1, t2)
#' }
elastic.partial.match <- function(template, target, argvals.template,
                                  argvals.target, lambda = 0,
                                  min.span = 0.3) {
  alignment_partial_match(as.numeric(template), as.numeric(target),
                          as.numeric(argvals.template),
                          as.numeric(argvals.target),
                          as.numeric(lambda), as.numeric(min.span))
}

# =============================================================================
# Curve Geodesic
# =============================================================================

#' Geodesic Path Between Two Curves
#'
#' Compute the geodesic (shortest path) between two curves in the elastic
#' metric. The geodesic is represented as a sequence of intermediate curves
#' interpolating between the two endpoints in the SRVF space.
#'
#' @param f1 Numeric vector of the first curve's values.
#' @param f2 Numeric vector of the second curve's values.
#' @param argvals Numeric vector of the evaluation grid (common to both curves).
#' @param n.points Number of points along the geodesic to compute (default 10).
#' @param lambda Regularisation parameter controlling warping smoothness
#'   (default 0).
#'
#' @return An object of class 'curve.geodesic' with components:
#' \describe{
#'   \item{path}{fdata of the interpolated curves along the geodesic}
#'   \item{distance}{the total geodesic (elastic) distance between f1 and f2}
#'   \item{gamma}{numeric vector of the optimal warping function}
#'   \item{argvals}{the evaluation grid}
#' }
#'
#' @references
#' Srivastava, A. and Klassen, E. (2016). \emph{Functional and Shape Data
#' Analysis}. Springer.
#'
#' @seealso \code{\link{elastic.distance}}
#'
#' @export
#' @examples
#' \donttest{
#' t <- seq(0, 1, length.out = 100)
#' f1 <- sin(2 * pi * t)
#' f2 <- cos(2 * pi * t)
#' geo <- curve.geodesic(f1, f2, t, n.points = 8)
#' }
curve.geodesic <- function(f1, f2, argvals, n.points = 10, lambda = 0) {
  res <- alignment_curve_geodesic(as.numeric(f1), as.numeric(f2),
                                  as.numeric(argvals),
                                  as.integer(n.points),
                                  as.numeric(lambda))

  result <- list(
    curves = fdata(res$curves, argvals = argvals),
    warps = res$warps,
    distances = res$distances,
    parameter_values = res$parameter_values,
    argvals = argvals
  )
  class(result) <- "curve.geodesic"
  result
}

#' @export
print.curve.geodesic <- function(x, ...) {
  n <- nrow(x$curves$data)
  cat("Curve Geodesic\n")
  cat("  Points along geodesic:", n, "\n")
  cat("  Total distance:", round(sum(x$distances), 4), "\n")
  invisible(x)
}

# =============================================================================
# Peak Persistence
# =============================================================================

#' Peak Persistence Analysis
#'
#' Analyse the persistence of peaks (local maxima) across a range of
#' regularisation strengths. As the regularisation parameter \eqn{\lambda}
#' increases, aligned curves become smoother and spurious peaks disappear.
#' This function tracks which peaks persist across the regularisation path,
#' providing a multi-scale view of the data's peak structure.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param lambdas Numeric vector of regularisation values to sweep. If NULL
#'   (default), uses \code{seq(0, 5, length.out = 30)}.
#' @param max.iter Maximum number of Karcher mean iterations at each lambda
#'   (default 10).
#' @param tol Convergence tolerance (default 1e-4).
#'
#' @return A list with components:
#' \describe{
#'   \item{lambdas}{the regularisation values used}
#'   \item{n_peaks}{integer vector of peak counts at each lambda}
#'   \item{peak_locations}{list of numeric vectors giving peak locations at each
#'     lambda}
#'   \item{persistence}{numeric vector of persistence scores for each peak}
#' }
#'
#' @references
#' Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
#' functional data using phase and amplitude separation.
#' \emph{Computational Statistics & Data Analysis}, 61:50--66.
#'
#' @seealso \code{\link{karcher.mean}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 100)
#' X <- matrix(0, 15, 100)
#' for (i in 1:15) X[i, ] <- sin(4 * pi * t + runif(1, -0.3, 0.3))
#' fd <- fdata(X, argvals = t)
#' pp <- peak.persistence(fd)
#' }
peak.persistence <- function(fdataobj, lambdas = NULL, max.iter = 10,
                             tol = 1e-4) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  if (is.null(lambdas)) {
    lambdas <- seq(0, 5, length.out = 30)
  }

  alignment_peak_persistence(fdataobj$data, fdataobj$argvals,
                             as.numeric(lambdas),
                             as.integer(max.iter),
                             as.numeric(tol))
}

# =============================================================================
# Transfer Alignment
# =============================================================================

#' Transfer Alignment Between Functional Samples
#'
#' Transfer the alignment learned from a source sample of curves to a target
#' sample. The Karcher mean and warping structure estimated from the source are
#' used to guide the alignment of the target curves, which is useful when the
#' target sample is small or noisy.
#'
#' @param source An object of class 'fdata' (the reference sample).
#' @param target An object of class 'fdata' (the sample to align).
#' @param lambda Regularisation parameter controlling warping smoothness
#'   (default 0).
#' @param max.iter Maximum number of iterations (default 15).
#' @param tol Convergence tolerance (default 1e-4).
#'
#' @return A list with components:
#' \describe{
#'   \item{aligned}{fdata of the aligned target curves}
#'   \item{gammas}{fdata of warping functions applied to the target}
#'   \item{source.mean}{fdata of the source Karcher mean}
#'   \item{distances}{numeric vector of elastic distances}
#' }
#'
#' @seealso \code{\link{karcher.mean}}, \code{\link{elastic.align}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 50)
#' X1 <- matrix(0, 20, 50)
#' for (i in 1:20) X1[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
#' X2 <- matrix(0, 5, 50)
#' for (i in 1:5) X2[i, ] <- sin(2 * pi * t + runif(1, -0.5, 0.5))
#' source <- fdata(X1, argvals = t)
#' target <- fdata(X2, argvals = t)
#' res <- transfer.alignment(source, target)
#' }
transfer.alignment <- function(source, target, lambda = 0, max.iter = 15,
                               tol = 1e-4) {
  if (!inherits(source, "fdata")) {
    stop("source must be of class 'fdata'")
  }
  if (!inherits(target, "fdata")) {
    stop("target must be of class 'fdata'")
  }

  argvals <- source$argvals
  res <- alignment_transfer(source$data, target$data, argvals,
                            as.numeric(lambda),
                            as.integer(max.iter),
                            as.numeric(tol))

  list(
    aligned = fdata(res$aligned_data, argvals = argvals),
    gammas = fdata(res$gammas, argvals = argvals),
    source.mean = fdata(matrix(res$source_mean, nrow = 1), argvals = argvals),
    distances = res$distances
  )
}

# =============================================================================
# Gaussian Generative Model
# =============================================================================

#' Gaussian Generative Model for Elastic Functional Data
#'
#' Fit a Gaussian generative model to elastic functional data using PCA on the
#' SRVF representations. The model captures amplitude variability and can
#' generate new random curves from the estimated distribution.
#'
#' @param karcher An object of class 'karcher.mean' (result of
#'   \code{\link{karcher.mean}}).
#' @param n.components Number of principal components to retain (default 3).
#' @param n.samples Number of random curves to generate (default 50).
#' @param seed Random seed for reproducibility (default 42).
#'
#' @return A list with components:
#' \describe{
#'   \item{samples}{fdata of generated random curves}
#'   \item{eigenvalues}{numeric vector of eigenvalues}
#'   \item{eigenfunctions}{matrix of eigenfunctions (components x grid)}
#'   \item{scores}{matrix of PCA scores (curves x components)}
#' }
#'
#' @references
#' Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
#' functional data using phase and amplitude separation.
#' \emph{Computational Statistics & Data Analysis}, 61:50--66.
#'
#' @seealso \code{\link{karcher.mean}}, \code{\link{joint.gauss.model}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 20, 50)
#' for (i in 1:20) X[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
#' fd <- fdata(X, argvals = t)
#' km <- karcher.mean(fd)
#' gm <- gauss.model(km, n.components = 2, n.samples = 10)
#' }
gauss.model <- function(karcher, n.components = 3, n.samples = 50,
                        seed = 42) {
  if (!inherits(karcher, "karcher.mean")) {
    stop("karcher must be of class 'karcher.mean'")
  }

  argvals <- karcher$mean$argvals
  res <- alignment_gauss_model(
    as.numeric(karcher$mean$data[1, ]),
    as.numeric(karcher$mean_srsf),
    karcher$gammas$data,
    karcher$aligned$data,
    argvals,
    as.integer(n.components),
    as.integer(n.samples),
    as.integer(seed)
  )

  list(
    samples = fdata(res$samples, argvals = argvals),
    eigenvalues = res$eigenvalues,
    eigenfunctions = res$eigenfunctions,
    scores = res$scores
  )
}

# =============================================================================
# Joint Gaussian Generative Model
# =============================================================================

#' Joint Amplitude-Phase Gaussian Generative Model
#'
#' Fit a joint Gaussian generative model that captures both amplitude and phase
#' variability. Unlike \code{\link{gauss.model}} which models amplitude only,
#' this model jointly models the amplitude and warping components, allowing
#' generated curves to exhibit realistic phase variability.
#'
#' @param karcher An object of class 'karcher.mean' (result of
#'   \code{\link{karcher.mean}}).
#' @param n.components Number of principal components to retain (default 3).
#' @param n.samples Number of random curves to generate (default 50).
#' @param balance Balance parameter controlling the relative weight of amplitude
#'   and phase components (default 1.0).
#' @param seed Random seed for reproducibility (default 42).
#'
#' @return A list with components:
#' \describe{
#'   \item{samples}{fdata of generated random curves}
#'   \item{eigenvalues}{numeric vector of eigenvalues}
#'   \item{eigenfunctions}{matrix of eigenfunctions}
#'   \item{scores}{matrix of joint PCA scores}
#' }
#'
#' @references
#' Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
#' functional data using phase and amplitude separation.
#' \emph{Computational Statistics & Data Analysis}, 61:50--66.
#'
#' @seealso \code{\link{karcher.mean}}, \code{\link{gauss.model}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 20, 50)
#' for (i in 1:20) X[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
#' fd <- fdata(X, argvals = t)
#' km <- karcher.mean(fd)
#' jgm <- joint.gauss.model(km, n.components = 2, n.samples = 10)
#' }
joint.gauss.model <- function(karcher, n.components = 3, n.samples = 50,
                              balance = 1.0, seed = 42) {
  if (!inherits(karcher, "karcher.mean")) {
    stop("karcher must be of class 'karcher.mean'")
  }

  argvals <- karcher$mean$argvals
  res <- alignment_joint_gauss_model(
    as.numeric(karcher$mean$data[1, ]),
    as.numeric(karcher$mean_srsf),
    karcher$gammas$data,
    karcher$aligned$data,
    argvals,
    as.integer(n.components),
    as.integer(n.samples),
    as.numeric(balance),
    as.integer(seed)
  )

  list(
    samples = fdata(res$samples, argvals = argvals),
    eigenvalues = res$eigenvalues,
    eigenfunctions = res$eigenfunctions,
    scores = res$scores
  )
}

# =============================================================================
# Horizontal fPCA (fPNS)
# =============================================================================

#' Horizontal Functional Principal Nested Spheres
#'
#' Perform horizontal functional PCA using the principal nested spheres (fPNS)
#' framework on warping functions. This captures the main modes of phase
#' variability from the warping functions estimated during Karcher mean
#' computation.
#'
#' @param karcher An object of class 'karcher.mean' (result of
#'   \code{\link{karcher.mean}}).
#' @param n.components Number of principal components to retain (default 3).
#'
#' @return A list with components:
#' \describe{
#'   \item{scores}{matrix of phase component scores (curves x components)}
#'   \item{eigenvalues}{numeric vector of eigenvalues}
#'   \item{eigenfunctions}{matrix of phase eigenfunctions (components x grid)}
#'   \item{variance.explained}{numeric vector of proportion of variance
#'     explained by each component}
#' }
#'
#' @references
#' Jung, S., Dryden, I.L., and Marron, J.S. (2012). Analysis of principal
#' nested spheres. \emph{Biometrika}, 99(3):551--568.
#'
#' Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
#' functional data using phase and amplitude separation.
#' \emph{Computational Statistics & Data Analysis}, 61:50--66.
#'
#' @seealso \code{\link{karcher.mean}}, \code{\link{gauss.model}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 20, 50)
#' for (i in 1:20) X[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
#' fd <- fdata(X, argvals = t)
#' km <- karcher.mean(fd)
#' hfp <- horiz.fpns(km, n.components = 2)
#' }
horiz.fpns <- function(karcher, n.components = 3) {
  if (!inherits(karcher, "karcher.mean")) {
    stop("karcher must be of class 'karcher.mean'")
  }

  argvals <- karcher$mean$argvals
  res <- alignment_horiz_fpns(
    as.numeric(karcher$mean$data[1, ]),
    as.numeric(karcher$mean_srsf),
    karcher$gammas$data,
    karcher$aligned$data,
    argvals,
    as.integer(n.components)
  )

  res
}
