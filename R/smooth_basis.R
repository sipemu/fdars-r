#' Penalized Basis Smoothing
#'
#' Smooth functional data using penalized B-spline or Fourier basis expansion.

#' Penalized Basis Smoothing
#'
#' Fits a penalized basis expansion to smooth functional data. Supports
#' B-spline and Fourier basis types with a roughness penalty.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param type Basis type: "bspline" or "fourier".
#' @param nbasis Number of basis functions. If NULL, defaults to min(ncol/4, 30).
#' @param lambda Smoothing parameter (penalty weight). Use 0 for no penalty.
#' @param lfd.order Order of the linear differential operator for the penalty
#'   (default 2 = penalize curvature).
#' @param period Period for Fourier basis (only used when type = "fourier").
#'   If NULL, defaults to the range of argvals.
#'
#' @return An object of class 'smooth.basis' with components:
#'   \item{coefficients}{Matrix of basis coefficients (nbasis x n)}
#'   \item{fitted}{Fitted smoothed curves as an fdata object}
#'   \item{edf}{Effective degrees of freedom}
#'   \item{gcv}{Generalized cross-validation score}
#'   \item{aic}{Akaike information criterion}
#'   \item{bic}{Bayesian information criterion}
#'   \item{nbasis}{Number of basis functions used}
#'   \item{lambda}{Smoothing parameter used}
#'   \item{type}{Basis type used}
#'   \item{fdataobj}{Original functional data}
#'
#' @export
smooth.basis.fd <- function(fdataobj, type = c("bspline", "fourier"),
                            nbasis = NULL, lambda = 0, lfd.order = 2,
                            period = NULL) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  type <- match.arg(type)
  m <- ncol(fdataobj$data)

  if (is.null(nbasis)) {
    nbasis <- min(max(m %/% 4, 5), 30)
  }

  if (is.null(period)) {
    period <- diff(range(fdataobj$argvals))
  }

  result <- .Call("wrap__smooth_basis_rust",
                  fdataobj$data, as.numeric(fdataobj$argvals),
                  type, as.integer(nbasis), as.numeric(lambda),
                  as.integer(lfd.order), as.numeric(period))

  if (is.null(result)) {
    stop("smooth.basis.fd failed: check data dimensions and parameters")
  }

  fitted_fd <- fdata(result$fitted, argvals = fdataobj$argvals,
                     names = fdataobj$names)

  structure(
    list(
      coefficients = result$coefficients,
      fitted = fitted_fd,
      edf = result$edf,
      gcv = result$gcv,
      aic = result$aic,
      bic = result$bic,
      nbasis = result$nbasis,
      lambda = lambda,
      type = type,
      fdataobj = fdataobj,
      call = match.call()
    ),
    class = "smooth.basis"
  )
}

#' Penalized Basis Smoothing with GCV-Optimal Lambda
#'
#' Selects the smoothing parameter by minimizing the generalized
#' cross-validation (GCV) criterion over a grid of lambda values.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param type Basis type: "bspline" or "fourier".
#' @param nbasis Number of basis functions. If NULL, defaults to min(ncol/4, 30).
#' @param lfd.order Order of the linear differential operator for the penalty.
#' @param log.lambda.range Range of log10(lambda) to search over.
#' @param n.grid Number of grid points for the lambda search.
#' @param period Period for Fourier basis.
#'
#' @return An object of class 'smooth.basis' (same as smooth.basis.fd).
#'
#' @export
smooth.basis.gcv <- function(fdataobj, type = c("bspline", "fourier"),
                             nbasis = NULL, lfd.order = 2,
                             log.lambda.range = c(-10, 5), n.grid = 50,
                             period = NULL) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  type <- match.arg(type)
  m <- ncol(fdataobj$data)

  if (is.null(nbasis)) {
    nbasis <- min(max(m %/% 4, 5), 30)
  }

  if (is.null(period)) {
    period <- diff(range(fdataobj$argvals))
  }

  result <- .Call("wrap__smooth_basis_gcv_rust",
                  fdataobj$data, as.numeric(fdataobj$argvals),
                  type, as.integer(nbasis), as.integer(lfd.order),
                  as.numeric(log.lambda.range[1]),
                  as.numeric(log.lambda.range[2]),
                  as.integer(n.grid), as.numeric(period))

  if (is.null(result)) {
    stop("smooth.basis.gcv failed: check data dimensions and parameters")
  }

  fitted_fd <- fdata(result$fitted, argvals = fdataobj$argvals,
                     names = fdataobj$names)

  structure(
    list(
      coefficients = result$coefficients,
      fitted = fitted_fd,
      edf = result$edf,
      gcv = result$gcv,
      aic = result$aic,
      bic = result$bic,
      nbasis = result$nbasis,
      lambda = 10^(result$gcv),  # Approximate from GCV-optimal
      type = type,
      fdataobj = fdataobj,
      call = match.call()
    ),
    class = "smooth.basis"
  )
}

#' @export
print.smooth.basis <- function(x, ...) {
  cat("Penalized Basis Smoothing\n")
  cat("  Basis type:", x$type, "\n")
  cat("  Number of basis functions:", x$nbasis, "\n")
  cat("  Lambda:", format(x$lambda, digits = 4), "\n")
  cat("  EDF:", format(x$edf, digits = 2), "\n")
  cat("  GCV:", format(x$gcv, digits = 6), "\n")
  invisible(x)
}

#' @export
plot.smooth.basis <- function(x, type = c("fit", "residuals", "gcv"), ...) {
  type <- match.arg(type)

  if (type == "fit") {
    n <- nrow(x$fdataobj$data)
    argvals <- x$fdataobj$argvals
    cols <- grDevices::adjustcolor(rep("grey60", n), alpha.f = 0.3)

    plot(NULL, xlim = range(argvals), ylim = range(x$fdataobj$data),
         xlab = x$fdataobj$names$xlab %||% "t",
         ylab = x$fdataobj$names$ylab %||% "x(t)",
         main = "Smooth Basis Fit", ...)
    for (i in seq_len(n)) {
      graphics::lines(argvals, x$fdataobj$data[i, ], col = cols[i])
    }
    for (i in seq_len(n)) {
      graphics::lines(argvals, x$fitted$data[i, ], col = "steelblue", lwd = 0.5)
    }
    graphics::legend("topright", legend = c("Raw", "Smoothed"),
                     col = c("grey60", "steelblue"), lty = 1, lwd = c(1, 1.5))
  } else if (type == "residuals") {
    resid <- x$fdataobj$data - x$fitted$data
    argvals <- x$fdataobj$argvals
    plot(NULL, xlim = range(argvals), ylim = range(resid),
         xlab = x$fdataobj$names$xlab %||% "t",
         ylab = "Residual",
         main = "Smoothing Residuals", ...)
    for (i in seq_len(nrow(resid))) {
      graphics::lines(argvals, resid[i, ],
                      col = grDevices::adjustcolor("steelblue", alpha.f = 0.3))
    }
    graphics::abline(h = 0, lty = 2)
  }

  invisible(x)
}
