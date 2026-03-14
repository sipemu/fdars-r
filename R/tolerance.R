#' Tolerance Bands for Functional Data
#'
#' Compute tolerance bands that are expected to contain a given fraction of
#' individual curves in the population.

#' Functional Tolerance Band
#'
#' Compute a tolerance band for functional data using one of several methods.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param method Method to use. One of "fpca" (default), "conformal", "scb",
#'   "exponential", "elastic", "phase", or "elastic.config".
#' @param coverage Target coverage probability (default 0.95).
#' @param ncomp Number of FPCA components (default 3). Used by "fpca",
#'   "exponential", "elastic", "phase", and "elastic.config" (as amplitude
#'   components for "elastic.config").
#' @param nb Number of bootstrap replicates (default 500). Used by "fpca",
#'   "scb", "exponential", "elastic", "phase", and "elastic.config".
#' @param band.type "pointwise" (default) or "simultaneous". Used by "fpca",
#'   "elastic", "phase", and "elastic.config".
#' @param cal.fraction Calibration fraction for conformal method (default 0.2).
#' @param score.type Nonconformity score: "supnorm" (default) or "l2".
#'   Used by "conformal".
#' @param bandwidth Kernel bandwidth for SCB Degras method. If NULL, a
#'   default is computed.
#' @param confidence Confidence level for SCB method (default is \code{coverage}).
#' @param multiplier Multiplier distribution: "gaussian" (default) or
#'   "rademacher". Used by "scb".
#' @param family Exponential family: "gaussian" (default), "binomial", or
#'   "poisson". Used by "exponential".
#' @param max.iter Maximum iterations for elastic method Karcher mean (default 10).
#' @param ncomp.phase Number of FPCA components for the phase band (default 3).
#'   Used by "elastic.config".
#' @param tol Convergence tolerance for Karcher mean (default 1e-4). Used by
#'   "elastic.config".
#' @param seed Random seed for reproducibility (default NULL).
#'
#' @return An object of class 'tolerance.band' with components:
#' \describe{
#'   \item{lower}{numeric vector of lower bounds}
#'   \item{upper}{numeric vector of upper bounds}
#'   \item{center}{numeric vector of center function}
#'   \item{half_width}{numeric vector of half-widths}
#'   \item{method}{the method used}
#'   \item{coverage}{the target coverage}
#'   \item{argvals}{evaluation points}
#'   \item{fdataobj}{the original fdata input}
#' }
#' The "phase" method additionally returns \code{gamma.lower},
#' \code{gamma.upper}, and \code{gamma.center} (warping function bounds),
#' with \code{lower}/\code{upper}/\code{center}/\code{half_width} from the
#' tangent-space band.
#'
#' The "elastic.config" method additionally returns \code{phase.lower},
#' \code{phase.upper}, and \code{phase.center} (warping function bounds),
#' with the primary \code{lower}/\code{upper}/\code{center}/\code{half_width}
#' from the amplitude band.
#'
#' Returns NULL with a warning if computation fails.
#'
#' @details
#' Available methods:
#' \describe{
#'   \item{fpca}{FPCA + bootstrap tolerance band. Reconstructs curves from PC
#'     scores and uses bootstrap to estimate pointwise or simultaneous quantiles.}
#'   \item{conformal}{Distribution-free conformal prediction band. Splits data
#'     into training and calibration sets.}
#'   \item{scb}{Simultaneous confidence band for the mean (Degras method).
#'     Uses multiplier bootstrap for critical values.}
#'   \item{exponential}{Tolerance band for exponential family functional data.
#'     Applies link function transformation.}
#'   \item{elastic}{Tolerance band in elastic (aligned) space. First computes
#'     Karcher mean, then applies FPCA band on aligned data.}
#'   \item{phase}{Phase tolerance band for warping function variation.
#'     Computes Karcher mean, extracts warping functions, and builds a
#'     tolerance band in the tangent (shooting vector) space. Returns both
#'     the warping-function bounds and the tangent-space band.}
#'   \item{elastic.config}{Joint amplitude and phase tolerance band with
#'     full configuration control. Separately controls the number of FPCA
#'     components for amplitude (\code{ncomp}) and phase (\code{ncomp.phase}),
#'     plus convergence tolerance (\code{tol}).}
#' }
#'
#' @references
#' Rathnayake, L.N. and Cuevas, A. (2016). Tolerance bands for functional data.
#' \emph{Technometrics}, 58(3):326--334.
#'
#' Lei, J. and Wasserman, L. (2014). Distribution-free prediction bands for
#' non-parametric regression. \emph{Journal of the Royal Statistical Society:
#' Series B}, 76(1):71--96.
#'
#' Degras, D. (2011). Simultaneous confidence bands for nonparametric
#' regression with functional data. \emph{Statistica Sinica}, 21(4):1735--1765.
#'
#' @export
#' @examples
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' \donttest{
#' band <- tolerance.band(fd, method = "fpca", nb = 100)
#' }
tolerance.band <- function(fdataobj,
                           method = c("fpca", "conformal", "scb",
                                      "exponential", "elastic",
                                      "phase", "elastic.config"),
                           coverage = 0.95,
                           ncomp = 3,
                           nb = 500,
                           band.type = c("pointwise", "simultaneous"),
                           cal.fraction = 0.2,
                           score.type = c("supnorm", "l2"),
                           bandwidth = NULL,
                           confidence = NULL,
                           multiplier = c("gaussian", "rademacher"),
                           family = c("gaussian", "binomial", "poisson"),
                           max.iter = 10,
                           ncomp.phase = 3,
                           tol = 1e-4,
                           seed = NULL) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  method <- match.arg(method)
  band.type <- match.arg(band.type)
  score.type <- match.arg(score.type)
  multiplier <- match.arg(multiplier)
  family <- match.arg(family)

  if (is.null(confidence)) confidence <- coverage

  data <- fdataobj$data
  argvals <- fdataobj$argvals
  seed_r <- if (is.null(seed)) NULL else as.integer(seed)

  res <- switch(method,
    "fpca" = tolerance_fpca(data, as.integer(ncomp), as.integer(nb),
                            coverage, band.type, seed_r),
    "conformal" = tolerance_conformal(data, cal.fraction, coverage,
                                      score.type, seed_r),
    "scb" = {
      bw <- if (is.null(bandwidth)) {
        diff(range(argvals)) / (length(argvals)^(1/5))
      } else {
        bandwidth
      }
      tolerance_scb_degras(data, argvals, bw, as.integer(nb),
                           confidence, multiplier)
    },
    "exponential" = tolerance_exponential(data, family, as.integer(ncomp),
                                          as.integer(nb), coverage, seed_r),
    "elastic" = tolerance_elastic(data, argvals, as.integer(ncomp),
                                  as.integer(nb), coverage, band.type,
                                  as.integer(max.iter), seed_r),
    "phase" = {
      res <- tolerance_phase_rust(data, argvals, as.integer(ncomp),
                                   as.integer(nb), coverage, band.type,
                                   as.integer(max.iter), seed_r)
      if (!res$success) {
        warning("Phase tolerance band computation failed")
        return(NULL)
      }
      result <- list(
        gamma.lower = res$gamma_lower,
        gamma.upper = res$gamma_upper,
        gamma.center = res$gamma_center,
        lower = res$tangent_lower,
        upper = res$tangent_upper,
        center = res$tangent_center,
        half_width = res$tangent_half_width,
        method = "phase",
        coverage = coverage,
        argvals = argvals,
        fdataobj = fdataobj
      )
      class(result) <- "tolerance.band"
      return(result)
    },
    "elastic.config" = {
      res <- tolerance_elastic_config_rust(data, argvals,
                                            as.integer(ncomp),
                                            as.integer(ncomp.phase),
                                            as.integer(nb), coverage,
                                            band.type,
                                            as.integer(max.iter), tol,
                                            seed_r)
      if (!res$success) {
        warning("Elastic config tolerance band computation failed")
        return(NULL)
      }
      result <- list(
        lower = res$amp_lower,
        upper = res$amp_upper,
        center = res$amp_center,
        half_width = res$amp_half_width,
        phase.lower = res$phase_gamma_lower,
        phase.upper = res$phase_gamma_upper,
        phase.center = res$phase_gamma_center,
        method = "elastic.config",
        coverage = coverage,
        argvals = argvals,
        fdataobj = fdataobj
      )
      class(result) <- "tolerance.band"
      return(result)
    }
  )

  if (!res$success) {
    warning("Tolerance band computation failed for method '", method, "'")
    return(NULL)
  }

  result <- list(
    lower = res$lower,
    upper = res$upper,
    center = res$center,
    half_width = res$half_width,
    method = method,
    coverage = coverage,
    argvals = argvals,
    fdataobj = fdataobj
  )
  class(result) <- "tolerance.band"
  result
}

# =============================================================================
# S3 methods
# =============================================================================

#' @export
print.tolerance.band <- function(x, ...) {
  cat("Functional Tolerance Band\n")
  cat("  Method:", x$method, "\n")
  cat("  Coverage:", x$coverage, "\n")
  cat("  Grid points:", length(x$argvals), "\n")
  cat("  Mean half-width:", round(mean(x$half_width), 4), "\n")
  invisible(x)
}

#' Plot Tolerance Band
#'
#' @param x An object of class 'tolerance.band'.
#' @param show.data Logical; if TRUE (default), overlay the original curves.
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object (invisibly).
#'
#' @export
plot.tolerance.band <- function(x, show.data = TRUE, ...) {
  m <- length(x$argvals)

  df_band <- data.frame(
    argval = x$argvals,
    lower = x$lower,
    upper = x$upper,
    center = x$center
  )

  p <- ggplot2::ggplot(df_band, ggplot2::aes(x = .data$argval))

  if (show.data && !is.null(x$fdataobj)) {
    n <- nrow(x$fdataobj$data)
    df_curves <- data.frame(
      curve_id = rep(seq_len(n), each = m),
      argval = rep(x$argvals, n),
      value = as.vector(t(x$fdataobj$data))
    )
    p <- p + ggplot2::geom_line(
      data = df_curves,
      ggplot2::aes(x = .data$argval, y = .data$value, group = .data$curve_id),
      alpha = 0.2, color = "grey60"
    )
  }

  p <- p +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      fill = "steelblue", alpha = 0.3
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$center),
      color = "steelblue", linewidth = 1
    ) +
    ggplot2::labs(
      title = paste0("Tolerance Band (", x$method, ", ",
                     x$coverage * 100, "% coverage)"),
      x = "t", y = "value"
    ) +
    ggplot2::theme_minimal()

  p
}
