#' Andrews Transformation
#'
#' Transforms a multivariate dataset into functional data using the Andrews
#' Fourier expansion. Each row of the input matrix becomes a curve on
#' \eqn{[-\pi, \pi]}{[-pi, pi]} defined by
#' \deqn{f_i(t) = x_{i1}/\sqrt{2} + x_{i2}\sin(t) + x_{i3}\cos(t) +
#'   x_{i4}\sin(2t) + x_{i5}\cos(2t) + \ldots}
#'
#' The resulting \code{fdata} object carries the Fourier basis matrix and
#' variable names as attributes, enabling [andrews_loadings()] to project
#' FPCA eigenfunctions back to the original variable space.
#'
#' @param X Numeric matrix or data.frame (\eqn{n \times p}). Rows are
#'   observations, columns are variables. Coerced via [as.matrix()].
#' @param m Number of grid points on \eqn{[-\pi, \pi]} (default 200).
#'
#' @return An [fdata] object (\eqn{n} curves, \eqn{m} grid points) with
#'   additional attributes:
#'   \describe{
#'     \item{\code{andrews_basis}}{The \eqn{m \times p} Fourier basis matrix.}
#'     \item{\code{andrews_varnames}}{Variable names (\code{colnames(X)} or
#'       \code{paste0("V", 1:p)} if \code{NULL}).}
#'   }
#'
#' @seealso [andrews_loadings()], [fdata2pc()]
#'
#' @export
#' @examples
#' X <- scale(iris[, 1:4])
#' fd <- andrews_transform(X)
#' plot(fd)
andrews_transform <- function(X, m = 200) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix or data.frame")
  }

  n <- nrow(X)
  p <- ncol(X)
  varnames <- colnames(X) %||% paste0("V", seq_len(p))
  t_grid <- seq(-pi, pi, length.out = m)

  # Build Fourier basis matrix [m x p]
  basis <- matrix(0, m, p)
  basis[, 1] <- 1 / sqrt(2)
  for (j in 2:p) {
    k <- ceiling((j - 1) / 2)
    if (j %% 2 == 0) {
      basis[, j] <- sin(k * t_grid)
    } else {
      basis[, j] <- cos(k * t_grid)
    }
  }

  fd <- fdata(X %*% t(basis), argvals = t_grid)
  attr(fd, "andrews_basis") <- basis
  attr(fd, "andrews_varnames") <- varnames
  fd
}


#' Andrews Loadings: Project FPCA Eigenfunctions to Original Variables
#'
#' Given FPCA results on Andrews-transformed data, projects each
#' eigenfunction back onto the Fourier basis to obtain loadings in the
#' original variable space. This reveals which variables contribute most
#' to each principal component.
#'
#' @param fpca Result from [fdata2pc()] (class \code{"fdata2pc"}).
#' @param fd_andrews The [fdata] object returned by [andrews_transform()]
#'   (must carry the \code{"andrews_basis"} attribute).
#' @param ncomp Number of principal components to project (default: all
#'   available in \code{fpca}).
#'
#' @return A data.frame of class \code{"andrews_loadings"} with columns:
#'   \describe{
#'     \item{\code{Variable}}{Original variable name.}
#'     \item{\code{Loading}}{Inner-product loading value.}
#'     \item{\code{PC}}{Principal component label (\code{"PC1"}, \code{"PC2"}, \ldots).}
#'   }
#'
#' @seealso [andrews_transform()], [fdata2pc()]
#'
#' @export
#' @examples
#' X <- scale(iris[, 1:4])
#' fd <- andrews_transform(X)
#' fpca <- fdata2pc(fd, ncomp = 3)
#' loadings <- andrews_loadings(fpca, fd)
#' head(loadings)
andrews_loadings <- function(fpca, fd_andrews, ncomp = NULL) {
  if (!inherits(fpca, "fdata2pc")) {
    stop("fpca must be a result from fdata2pc() (class 'fdata2pc')")
  }

  basis <- attr(fd_andrews, "andrews_basis")
  if (is.null(basis)) {
    stop("fd_andrews must be created by andrews_transform() ",
         "(missing 'andrews_basis' attribute)")
  }

  varnames <- attr(fd_andrews, "andrews_varnames")
  n_available <- nrow(fpca$rotation$data)
  if (is.null(ncomp)) ncomp <- n_available
  ncomp <- min(ncomp, n_available)

  t_grid <- fd_andrews$argvals
  dt <- diff(t_grid[1:2])

  loading_list <- vector("list", ncomp)
  for (pc_idx in seq_len(ncomp)) {
    eigen_vals <- fpca$rotation$data[pc_idx, ]
    loadings <- as.numeric(t(basis) %*% eigen_vals) * dt
    loading_list[[pc_idx]] <- data.frame(
      Variable = varnames,
      Loading = loadings,
      PC = paste0("PC", pc_idx),
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, loading_list)
  rownames(result) <- NULL
  class(result) <- c("andrews_loadings", "data.frame")
  result
}
