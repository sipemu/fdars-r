#' Scalar-on-Shape Regression
#'
#' Fits a scalar-on-shape regression model where the response depends on the
#' \emph{shape} of the functional predictor, with phase variation removed via
#' elastic alignment (Fisher-Rao framework).
#'
#' @param fdataobj An object of class 'fdata'.
#' @param y Scalar response vector (length n).
#' @param nbasis Number of Fourier basis functions for beta (default 11).
#' @param lambda Roughness penalty weight (default 1e-3).
#' @param index.method Index function method: "identity" (default),
#'   "polynomial", or "nadaraya.watson".
#' @param index.param Parameter for the index method: polynomial degree
#'   (integer) or Nadaraya-Watson bandwidth (numeric). Default 2.
#' @param g.degree Polynomial degree for the amplitude link g (default 2).
#' @param dp.lambda Dynamic programming alignment penalty (default 0).
#' @param max.iter.outer Maximum outer iterations (default 10).
#' @param max.iter.inner Maximum inner iterations (default 15).
#' @param tol Convergence tolerance (default 1e-4).
#'
#' @return An object of class 'scalar.on.shape' (a list) with components:
#'   \item{beta}{Estimated shape index function (length m).}
#'   \item{beta.coefficients}{Fourier coefficients for beta.}
#'   \item{gammas}{Warping functions (n x m matrix).}
#'   \item{shape.scores}{Shape scores (length n).}
#'   \item{h.coefficients}{Index link function coefficients.}
#'   \item{g.coefficients}{Amplitude link function coefficients.}
#'   \item{fitted.values}{Fitted response values.}
#'   \item{residuals}{Residuals.}
#'   \item{sse}{Residual sum of squares.}
#'   \item{r.squared}{Coefficient of determination.}
#'   \item{n.iter.outer}{Number of outer iterations.}
#'   \item{n.iter.inner}{Number of inner iterations.}
#'   \item{index.method}{Index method used.}
#'   \item{index.param}{Index method parameter.}
#'   \item{fdataobj}{Training fdata object.}
#'   \item{y}{Training response.}
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' fit <- scalar.on.shape(fd, y, nbasis = 5)
#' }
#'
#' @export
scalar.on.shape <- function(fdataobj, y, nbasis = 11, lambda = 1e-3,
                            index.method = c("identity", "polynomial",
                                             "nadaraya.watson"),
                            index.param = 2, g.degree = 2, dp.lambda = 0.0,
                            max.iter.outer = 10, max.iter.inner = 15,
                            tol = 1e-4) {
  if (!inherits(fdataobj, "fdata")) stop("fdataobj must be of class 'fdata'")
  index.method <- match.arg(index.method)
  argvals <- as.numeric(fdataobj$argvals)

  result <- .Call("wrap__scalar_on_shape_rust",
                  fdataobj$data, as.numeric(y), argvals,
                  as.integer(nbasis), as.numeric(lambda),
                  index.method, as.numeric(index.param),
                  as.integer(g.degree), as.numeric(dp.lambda),
                  as.integer(max.iter.outer), as.integer(max.iter.inner),
                  as.numeric(tol))

  if (is.null(result)) stop("scalar.on.shape failed")

  out <- list(
    beta = result$beta,
    beta.coefficients = result$beta_coefficients,
    gammas = result$gammas,
    shape.scores = result$shape_scores,
    h.coefficients = result$h_coefficients,
    g.coefficients = result$g_coefficients,
    fitted.values = result$fitted_values,
    residuals = result$residuals,
    sse = result$sse,
    r.squared = result$r_squared,
    n.iter.outer = result$n_iter_outer,
    n.iter.inner = result$n_iter_inner,
    index.method = result$index_method,
    index.param = result$index_param,
    fdataobj = fdataobj,
    y = y
  )
  class(out) <- "scalar.on.shape"
  out
}

#' Predict from a Scalar-on-Shape Model
#'
#' @param object A fitted \code{scalar.on.shape} model.
#' @param newdata An object of class 'fdata'.
#' @param ... Ignored.
#'
#' @return Predicted response values (numeric vector).
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' fit <- scalar.on.shape(fd, y, nbasis = 5)
#' pred <- predict(fit, fd[1:10, ])
#' }
#'
#' @export
predict.scalar.on.shape <- function(object, newdata, ...) {
  if (!inherits(newdata, "fdata")) stop("newdata must be of class 'fdata'")
  argvals <- as.numeric(object$fdataobj$argvals)

  result <- .Call("wrap__predict_scalar_on_shape_rust",
                  as.numeric(object$beta),
                  as.numeric(object$beta.coefficients),
                  object$gammas,
                  as.numeric(object$shape.scores),
                  as.numeric(object$h.coefficients),
                  as.numeric(object$g.coefficients),
                  as.numeric(object$fitted.values),
                  as.numeric(object$residuals),
                  as.numeric(object$sse),
                  as.numeric(object$r.squared),
                  as.integer(object$n.iter.outer),
                  as.integer(object$n.iter.inner),
                  object$index.method,
                  as.numeric(object$index.param),
                  newdata$data, argvals)

  if (is.null(result)) stop("predict.scalar.on.shape failed")
  result
}

#' Print Scalar-on-Shape Model
#' @param x A \code{scalar.on.shape} object.
#' @param ... Ignored.
#' @export
print.scalar.on.shape <- function(x, ...) {
  cat("Scalar-on-Shape Regression\n")
  cat("  Index method:", x$index.method, "\n")
  cat("  R-squared:", round(x$r.squared, 4), "\n")
  cat("  SSE:", round(x$sse, 4), "\n")
  cat("  Iterations:", x$n.iter.outer, "(outer),", x$n.iter.inner, "(inner)\n")
  invisible(x)
}
