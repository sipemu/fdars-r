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
elastic.align <- function(fdataobj, target = NULL, periodic = FALSE) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  argvals <- fdataobj$argvals
  rotations <- NULL

  # Periodic pre-processing: circularly rotate curves to canonical position
  if (periodic) {
    rot <- .periodic_rotate(fdataobj$data)
    work_data <- rot$data
    rotations <- rot$shifts

    if (!is.null(target)) {
      if (inherits(target, "fdata")) {
        target <- as.numeric(target$data[1, ])
      }
      target_rot <- .periodic_rotate(matrix(target, nrow = 1))
      target <- as.numeric(target_rot$data[1, ])
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
    rotations = rotations
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
karcher.mean <- function(fdataobj, max.iter = 20, tol = 1e-4, periodic = FALSE) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  argvals <- fdataobj$argvals
  rotations <- NULL

  if (periodic) {
    rot <- .periodic_rotate(fdataobj$data)
    work_data <- rot$data
    rotations <- rot$shifts
  } else {
    work_data <- fdataobj$data
  }

  res <- alignment_karcher_mean(work_data, argvals,
                                as.integer(max.iter), as.numeric(tol))

  result <- list(
    mean = fdata(matrix(res$mean, nrow = 1), argvals = argvals),
    aligned = fdata(res$aligned_data, argvals = argvals),
    gammas = fdata(res$gammas, argvals = argvals),
    n.iter = res$n_iter,
    converged = res$converged,
    fdataobj = fdataobj,
    rotations = rotations
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
    cat("  Periodic: TRUE\n")
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
    cat("  Periodic: TRUE\n")
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
#' Circularly rotate each curve so that its global maximum is at a canonical
#' grid position. This is useful as a pre-processing step before elastic
#' alignment of periodic functional data (e.g., data on \eqn{[0, 2\pi]} where
#' \eqn{f(0) = f(2\pi)}).
#'
#' @param fdataobj An object of class 'fdata'.
#' @param target.pos Target grid index for the global maximum. If NULL
#'   (default), uses \code{floor(ncol(fdataobj$data) / 4)}.
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
periodic.rotate <- function(fdataobj, target.pos = NULL) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }
  rot <- .periodic_rotate(fdataobj$data, target_pos = target.pos)
  list(
    fdataobj = fdata(rot$data, argvals = fdataobj$argvals),
    shifts = rot$shifts
  )
}

#' @noRd
.periodic_rotate <- function(data_matrix, target_pos = NULL) {
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
    # Circular shift: move elements so that max_idx lands at target_pos
    idx <- ((seq_len(m) - 1L + shift) %% m) + 1L
    rotated[i, ] <- data_matrix[i, idx]
  }

  list(data = rotated, shifts = shifts)
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
