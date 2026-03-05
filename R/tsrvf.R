#' TSRVF: Transported Square-Root Velocity Function
#'
#' Functions for mapping aligned functional data to a tangent space
#' on the Hilbert sphere, enabling standard linear operations (PCA,
#' regression) on elastically aligned curves.

#' TSRVF Transform
#'
#' Compute the TSRVF (Transported SRSF) representation of functional data.
#' This performs Karcher mean alignment and then maps the aligned SRSFs
#' to the tangent space at the mean via the inverse exponential map.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param max.iter Maximum number of Karcher mean iterations. Default 20.
#' @param tol Convergence tolerance for Karcher mean. Default 1e-4.
#' @param lambda Penalty weight on warp deviation from identity. Default 0.
#'
#' @return An object of class 'tsrvf' with components:
#' \describe{
#'   \item{tangent_vectors}{fdata of tangent vectors in Euclidean space}
#'   \item{mean}{fdata of the Karcher mean curve}
#'   \item{gammas}{fdata of warping functions}
#'   \item{converged}{logical indicating Karcher mean convergence}
#'   \item{fdataobj}{the original input}
#' }
#'
#' @references
#' Su, J., Kurtek, S., Klassen, E., and Srivastava, A. (2014).
#' Statistical analysis of trajectories on Riemannian manifolds:
#' bird migration, hurricane tracking and video surveillance.
#' \emph{Annals of Applied Statistics}, 8(1):530--552.
#'
#' @export
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
#' tv <- tsrvf.transform(fd, max.iter = 5)
#' }
tsrvf.transform <- function(fdataobj, max.iter = 20, tol = 1e-4, lambda = 0) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  argvals <- fdataobj$argvals
  res <- alignment_tsrvf_transform(fdataobj$data, as.numeric(argvals),
                                   as.integer(max.iter), as.numeric(tol),
                                   as.numeric(lambda))

  result <- list(
    tangent_vectors = fdata(res$tangent_vectors, argvals = argvals),
    mean = fdata(matrix(res$mean, nrow = 1), argvals = argvals),
    mean_srsf = res$mean_srsf,
    gammas = fdata(res$gammas, argvals = argvals),
    converged = res$converged,
    fdataobj = fdataobj
  )
  class(result) <- "tsrvf"
  result
}

#' TSRVF from Pre-computed Alignment
#'
#' Compute the TSRVF representation from a pre-computed Karcher mean,
#' avoiding the expensive re-computation of the alignment.
#'
#' @param karcher An object of class 'karcher.mean'.
#'
#' @return An object of class 'tsrvf'.
#'
#' @export
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
#' km <- karcher.mean(fd, max.iter = 5)
#' tv <- tsrvf.from.alignment(km)
#' }
tsrvf.from.alignment <- function(karcher) {
  if (!inherits(karcher, "karcher.mean")) {
    stop("karcher must be of class 'karcher.mean'")
  }

  argvals <- karcher$mean$argvals
  mean_vec <- as.numeric(karcher$mean$data[1, ])
  mean_srsf <- as.numeric(srsf.transform(karcher$mean)$data[1, ])

  res <- alignment_tsrvf_from_karcher(
    mean_vec, mean_srsf,
    karcher$gammas$data, karcher$aligned$data,
    as.numeric(argvals),
    as.integer(karcher$n.iter), karcher$converged
  )

  result <- list(
    tangent_vectors = fdata(res$tangent_vectors, argvals = argvals),
    mean = karcher$mean,
    mean_srsf = res$mean_srsf,
    gammas = karcher$gammas,
    converged = res$converged,
    fdataobj = karcher$fdataobj
  )
  class(result) <- "tsrvf"
  result
}

#' Inverse TSRVF Transform
#'
#' Reconstruct aligned curves from their TSRVF tangent vector representation.
#'
#' @param tsrvf.obj An object of class 'tsrvf'.
#'
#' @return An fdata object of reconstructed aligned curves.
#'
#' @export
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
#' tv <- tsrvf.transform(fd, max.iter = 5)
#' recon <- tsrvf.inverse(tv)
#' }
tsrvf.inverse <- function(tsrvf.obj) {
  if (!inherits(tsrvf.obj, "tsrvf")) {
    stop("tsrvf.obj must be of class 'tsrvf'")
  }

  argvals <- tsrvf.obj$mean$argvals
  mean_vec <- as.numeric(tsrvf.obj$mean$data[1, ])
  mean_srsf <- tsrvf.obj$mean_srsf

  # Compute mean_srsf_norm and srsf_norms
  # These are L2 norms on [0,1]
  m <- length(argvals)
  time_01 <- seq(0, 1, length.out = m)
  h <- 1 / (m - 1)

  mean_srsf_norm <- sqrt(sum(mean_srsf^2) * h)

  n <- nrow(tsrvf.obj$tangent_vectors$data)
  # Approximate srsf_norms as mean_srsf_norm (close for small tangent vectors)
  srsf_norms <- rep(mean_srsf_norm, n)

  result <- alignment_tsrvf_inverse(
    tsrvf.obj$tangent_vectors$data, mean_vec,
    mean_srsf, mean_srsf_norm,
    srsf_norms, tsrvf.obj$gammas$data,
    as.numeric(argvals), tsrvf.obj$converged
  )

  fdata(result, argvals = argvals)
}

#' @export
print.tsrvf <- function(x, ...) {
  n <- nrow(x$tangent_vectors$data)
  m <- ncol(x$tangent_vectors$data)
  cat("TSRVF (Transported SRSF)\n")
  cat("  Curves:", n, "x", m, "grid points\n")
  cat("  Converged:", x$converged, "\n")
  invisible(x)
}

#' Plot TSRVF Results
#'
#' @param x An object of class 'tsrvf'.
#' @param type Character: "tangent" (default), "mean", or "gammas".
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object (invisibly).
#'
#' @export
plot.tsrvf <- function(x, type = c("tangent", "mean", "gammas"), ...) {
  type <- match.arg(type)

  if (type == "gammas") {
    p <- .plot_fdata_curves(x$gammas, title = "Warping Functions")
  } else if (type == "mean") {
    n <- nrow(x$tangent_vectors$data)
    m <- ncol(x$tangent_vectors$data)
    argvals <- x$mean$argvals

    df_mean <- data.frame(
      argval = argvals,
      value = as.numeric(x$mean$data[1, ])
    )

    p <- ggplot2::ggplot(df_mean, ggplot2::aes(x = .data$argval, y = .data$value)) +
      ggplot2::geom_line(color = "red", linewidth = 1.2) +
      ggplot2::labs(title = "Karcher Mean", x = "t", y = "value") +
      ggplot2::theme_minimal()
  } else {
    p <- .plot_fdata_curves(x$tangent_vectors, title = "TSRVF Tangent Vectors")
  }

  p
}
