#' Elastic FPCA
#'
#' Amplitude, phase, and joint FPCA after elastic alignment using the
#' Square Root Velocity Framework (SRVF).

#' Vertical (Amplitude) FPCA
#'
#' Performs FPCA on the amplitude component of curves after elastic alignment.
#' Captures shape variability independent of timing/phase.
#'
#' @param karcher A Karcher mean result from \code{\link{karcher.mean}}.
#' @param ncomp Number of principal components (default 3).
#'
#' @return An object of class 'vert.fpca' with components:
#'   \item{scores}{Matrix of FPC scores (n x ncomp)}
#'   \item{eigenfunctions.q}{SRVF eigenfunctions}
#'   \item{eigenfunctions.f}{Original-space eigenfunctions}
#'   \item{eigenvalues}{Eigenvalues (variance explained by each component)}
#'   \item{cumulative.variance}{Cumulative proportion of variance}
#'   \item{mean.q}{Mean SRVF}
#'   \item{karcher}{The Karcher mean used}
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' km <- karcher.mean(fd)
#' v <- vert.fpca(km, ncomp = 2)
#' v
#' }
#'
#' @export
vert.fpca <- function(karcher, ncomp = 3) {
  if (!inherits(karcher, "karcher.mean")) {
    stop("karcher must be of class 'karcher.mean'")
  }

  result <- .Call("wrap__elastic_vert_fpca_rust",
                  karcher$aligned$data, karcher$gammas$data,
                  as.numeric(karcher$mean$data[1, ]),
                  as.numeric(karcher$mean_srsf),
                  NULL,
                  as.numeric(karcher$mean$argvals),
                  as.integer(ncomp))

  if (is.null(result)) {
    stop("vert.fpca failed: check karcher mean result and ncomp")
  }

  structure(
    list(
      scores = result$scores,
      eigenfunctions.q = result$eigenfunctions_q,
      eigenfunctions.f = result$eigenfunctions_f,
      eigenvalues = result$eigenvalues,
      cumulative.variance = result$cumulative_variance,
      mean.q = result$mean_q,
      karcher = karcher,
      ncomp = ncomp,
      call = match.call()
    ),
    class = "vert.fpca"
  )
}

#' Horizontal (Phase) FPCA
#'
#' Performs FPCA on the phase (warping) component after elastic alignment.
#' Captures timing variability independent of shape.
#'
#' @param karcher A Karcher mean result from \code{\link{karcher.mean}}.
#' @param ncomp Number of principal components (default 3).
#'
#' @return An object of class 'horiz.fpca' with components:
#'   \item{scores}{Matrix of FPC scores (n x ncomp)}
#'   \item{eigenfunctions.psi}{Eigenfunctions in tangent space}
#'   \item{eigenfunctions.gam}{Eigenfunctions as warp functions}
#'   \item{eigenvalues}{Eigenvalues}
#'   \item{cumulative.variance}{Cumulative proportion of variance}
#'   \item{mean.psi}{Mean tangent vector}
#'   \item{shooting.vectors}{Shooting vectors matrix}
#'   \item{karcher}{The Karcher mean used}
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' km <- karcher.mean(fd)
#' h <- horiz.fpca(km, ncomp = 2)
#' h
#' }
#'
#' @export
horiz.fpca <- function(karcher, ncomp = 3) {
  if (!inherits(karcher, "karcher.mean")) {
    stop("karcher must be of class 'karcher.mean'")
  }

  result <- .Call("wrap__elastic_horiz_fpca_rust",
                  karcher$aligned$data, karcher$gammas$data,
                  as.numeric(karcher$mean$data[1, ]),
                  as.numeric(karcher$mean_srsf),
                  NULL,
                  as.numeric(karcher$mean$argvals),
                  as.integer(ncomp))

  if (is.null(result)) {
    stop("horiz.fpca failed: check karcher mean result and ncomp")
  }

  structure(
    list(
      scores = result$scores,
      eigenfunctions.psi = result$eigenfunctions_psi,
      eigenfunctions.gam = result$eigenfunctions_gam,
      eigenvalues = result$eigenvalues,
      cumulative.variance = result$cumulative_variance,
      mean.psi = result$mean_psi,
      shooting.vectors = result$shooting_vectors,
      karcher = karcher,
      ncomp = ncomp,
      call = match.call()
    ),
    class = "horiz.fpca"
  )
}

#' Joint (Amplitude + Phase) FPCA
#'
#' Performs joint FPCA combining amplitude and phase components with a
#' balancing parameter.
#'
#' @param karcher A Karcher mean result from \code{\link{karcher.mean}}.
#' @param ncomp Number of principal components (default 3).
#'
#' @return An object of class 'joint.fpca' with components:
#'   \item{scores}{Matrix of joint FPC scores (n x ncomp)}
#'   \item{eigenvalues}{Eigenvalues}
#'   \item{cumulative.variance}{Cumulative proportion of variance}
#'   \item{balance.c}{Balancing parameter between amplitude and phase}
#'   \item{vert.component}{Vertical (amplitude) eigenvector components}
#'   \item{horiz.component}{Horizontal (phase) eigenvector components}
#'   \item{karcher}{The Karcher mean used}
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' km <- karcher.mean(fd)
#' j <- joint.fpca(km, ncomp = 2)
#' j
#' }
#'
#' @export
joint.fpca <- function(karcher, ncomp = 3) {
  if (!inherits(karcher, "karcher.mean")) {
    stop("karcher must be of class 'karcher.mean'")
  }

  result <- .Call("wrap__elastic_joint_fpca_rust",
                  karcher$aligned$data, karcher$gammas$data,
                  as.numeric(karcher$mean$data[1, ]),
                  as.numeric(karcher$mean_srsf),
                  NULL,
                  as.numeric(karcher$mean$argvals),
                  as.integer(ncomp))

  if (is.null(result)) {
    stop("joint.fpca failed: check karcher mean result and ncomp")
  }

  structure(
    list(
      scores = result$scores,
      eigenvalues = result$eigenvalues,
      cumulative.variance = result$cumulative_variance,
      balance.c = result$balance_c,
      vert.component = result$vert_component,
      horiz.component = result$horiz_component,
      karcher = karcher,
      ncomp = ncomp,
      call = match.call()
    ),
    class = "joint.fpca"
  )
}

#' @export
print.vert.fpca <- function(x, ...) {
  cat("Vertical (Amplitude) FPCA\n")
  cat("  Components:", x$ncomp, "\n")
  cat("  Cumulative variance:", paste(format(x$cumulative.variance * 100, digits = 1), "%", collapse = ", "), "\n")
  invisible(x)
}

#' @export
print.horiz.fpca <- function(x, ...) {
  cat("Horizontal (Phase) FPCA\n")
  cat("  Components:", x$ncomp, "\n")
  cat("  Cumulative variance:", paste(format(x$cumulative.variance * 100, digits = 1), "%", collapse = ", "), "\n")
  invisible(x)
}

#' @export
print.joint.fpca <- function(x, ...) {
  cat("Joint (Amplitude + Phase) FPCA\n")
  cat("  Components:", x$ncomp, "\n")
  cat("  Balance parameter:", format(x$balance.c, digits = 3), "\n")
  cat("  Cumulative variance:", paste(format(x$cumulative.variance * 100, digits = 1), "%", collapse = ", "), "\n")
  invisible(x)
}

#' @export
plot.vert.fpca <- function(x, type = c("scores", "eigenfunctions", "variance"), ...) {
  type <- match.arg(type)
  if (type == "variance") {
    df <- data.frame(
      PC = factor(paste0("PC", seq_along(x$eigenvalues)),
                  levels = paste0("PC", seq_along(x$eigenvalues))),
      Variance = x$eigenvalues / sum(x$eigenvalues) * 100
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PC, y = .data$Variance)) +
      ggplot2::geom_col(fill = "#4A90D9", alpha = 0.7) +
      ggplot2::labs(title = "Vertical FPCA - Scree Plot", x = NULL, y = "% Variance") +
      ggplot2::theme_minimal()
    print(p)
  } else if (type == "scores") {
    if (ncol(x$scores) >= 2) {
      df <- data.frame(PC1 = x$scores[, 1], PC2 = x$scores[, 2])
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PC1, y = .data$PC2)) +
        ggplot2::geom_point(color = "#4A90D9", size = 2.5, alpha = 0.7) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
        ggplot2::labs(title = "Vertical FPCA Scores", x = "PC1", y = "PC2") +
        ggplot2::theme_minimal()
      print(p)
    }
  } else if (type == "eigenfunctions") {
    nc <- min(ncol(x$eigenfunctions.f), x$ncomp)
    nr <- nrow(x$eigenfunctions.f)
    argvals <- as.numeric(x$karcher$mean$argvals)
    if (length(argvals) != nr) argvals <- seq_len(nr)
    df <- data.frame(
      t = rep(argvals, nc),
      value = as.vector(x$eigenfunctions.f[, seq_len(nc), drop = FALSE]),
      PC = factor(rep(paste0("PC", seq_len(nc)), each = nr))
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$t, y = .data$value)) +
      ggplot2::geom_line(color = "#4A90D9", linewidth = 1) +
      ggplot2::facet_wrap(~ .data$PC, scales = "free_y") +
      ggplot2::labs(title = "Amplitude Eigenfunctions", x = "t", y = "Eigenfunction") +
      ggplot2::theme_minimal()
    print(p)
  }
  invisible(x)
}

#' @export
plot.horiz.fpca <- function(x, type = c("scores", "eigenfunctions", "warps", "variance"), ...) {
  type <- match.arg(type)
  if (type == "variance") {
    df <- data.frame(
      PC = factor(paste0("PC", seq_along(x$eigenvalues)),
                  levels = paste0("PC", seq_along(x$eigenvalues))),
      Variance = x$eigenvalues / sum(x$eigenvalues) * 100
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PC, y = .data$Variance)) +
      ggplot2::geom_col(fill = "#E6A020", alpha = 0.7) +
      ggplot2::labs(title = "Horizontal FPCA - Scree Plot", x = NULL, y = "% Variance") +
      ggplot2::theme_minimal()
    print(p)
  } else if (type == "scores") {
    if (ncol(x$scores) >= 2) {
      df <- data.frame(PC1 = x$scores[, 1], PC2 = x$scores[, 2])
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PC1, y = .data$PC2)) +
        ggplot2::geom_point(color = "#E6A020", size = 2.5, alpha = 0.7) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
        ggplot2::labs(title = "Phase FPCA Scores", x = "PC1", y = "PC2") +
        ggplot2::theme_minimal()
      print(p)
    }
  } else if (type == "warps" || type == "eigenfunctions") {
    nc <- min(ncol(x$eigenfunctions.gam), x$ncomp)
    nr <- nrow(x$eigenfunctions.gam)
    argvals <- as.numeric(x$karcher$mean$argvals)
    if (length(argvals) != nr) argvals <- seq_len(nr)
    df <- data.frame(
      t = rep(argvals, nc),
      value = as.vector(x$eigenfunctions.gam[, seq_len(nc), drop = FALSE]),
      PC = factor(rep(paste0("PC", seq_len(nc)), each = nr))
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$t, y = .data$value)) +
      ggplot2::geom_line(color = "#E6A020", linewidth = 1) +
      ggplot2::facet_wrap(~ .data$PC, scales = "free_y") +
      ggplot2::labs(title = "Phase Eigenfunctions", x = "t", y = "Warp Eigenfunction") +
      ggplot2::theme_minimal()
    print(p)
  }
  invisible(x)
}

#' @export
plot.joint.fpca <- function(x, type = c("scores", "variance", "balance"), ...) {
  type <- match.arg(type)
  if (type == "variance") {
    df <- data.frame(
      PC = factor(paste0("PC", seq_along(x$eigenvalues)),
                  levels = paste0("PC", seq_along(x$eigenvalues))),
      Variance = x$eigenvalues / sum(x$eigenvalues) * 100
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PC, y = .data$Variance)) +
      ggplot2::geom_col(fill = "#7B2D8E", alpha = 0.7) +
      ggplot2::labs(title = "Joint FPCA - Scree Plot", x = NULL, y = "% Variance") +
      ggplot2::theme_minimal()
    print(p)
  } else if (type == "scores") {
    if (ncol(x$scores) >= 2) {
      df <- data.frame(PC1 = x$scores[, 1], PC2 = x$scores[, 2])
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PC1, y = .data$PC2)) +
        ggplot2::geom_point(color = "#7B2D8E", size = 2.5, alpha = 0.7) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
        ggplot2::labs(title = "Joint FPCA Scores", x = "PC1", y = "PC2") +
        ggplot2::theme_minimal()
      print(p)
    }
  } else if (type == "balance") {
    npc <- ncol(x$vert.component)
    vert_var <- colSums(x$vert.component^2)
    horiz_var <- colSums(x$horiz.component^2)
    total <- vert_var + horiz_var
    df <- data.frame(
      PC = rep(factor(paste0("PC", seq_len(npc)),
                      levels = paste0("PC", seq_len(npc))), 2),
      Proportion = c(vert_var / total, horiz_var / total),
      Source = factor(rep(c("Amplitude", "Phase"), each = npc))
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PC, y = .data$Proportion,
                                           fill = .data$Source)) +
      ggplot2::geom_col(position = "stack", alpha = 0.7) +
      ggplot2::scale_fill_manual(values = c("Amplitude" = "#4A90D9", "Phase" = "#E6A020")) +
      ggplot2::labs(title = paste0("Joint FPCA Balance (c = ",
                                    format(x$balance.c, digits = 3), ")"),
                    x = NULL, y = "Proportion", fill = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")
    print(p)
  }
  invisible(x)
}
