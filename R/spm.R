#' Statistical Process Monitoring for Functional Data
#'
#' Functions for monitoring functional data processes using FPCA-based
#' control charts, supporting both univariate and multivariate settings.
#'
#' @name spm
NULL

# =============================================================================
# Phase I: Build control chart
# =============================================================================

#' Build Univariate SPM Control Chart (Phase I)
#'
#' Builds an FPCA-based control chart from in-control functional data.
#' The data is split into a tuning set (for FPCA) and a calibration set
#' (for control limits). Computes Hotelling T-squared and SPE statistics.
#'
#' @param fdataobj An object of class \code{fdata} containing in-control
#'   functional data (at least 4 observations).
#' @param ncomp Number of principal components (default 5).
#' @param alpha Significance level for control limits (default 0.05).
#' @param tuning.fraction Fraction of data used for FPCA tuning (default 0.5).
#' @param seed Random seed for train/calibration split (default 42).
#'
#' @return An object of class \code{spm.chart} with components:
#' \describe{
#'   \item{eigenvalues}{Eigenvalues from FPCA}
#'   \item{t2.phase1}{T-squared values for calibration set}
#'   \item{spe.phase1}{SPE values for calibration set}
#'   \item{t2.ucl}{T-squared upper control limit}
#'   \item{spe.ucl}{SPE upper control limit}
#'   \item{ncomp}{Number of components used}
#'   \item{t2.description}{Description of T-squared limit}
#'   \item{spe.description}{Description of SPE limit}
#'   \item{fdataobj}{Original fdata object}
#'   \item{.rust}{Internal fields for Phase II monitoring}
#' }
#'
#' @seealso \code{\link{spm.monitor}} for Phase II monitoring,
#'   \code{\link{spm.ewma}} for EWMA-based monitoring,
#'   \code{\link{frcc.phase1}} for regression-adjusted monitoring
#'
#' @examples
#' \donttest{
#' # Simulate in-control functional data
#' set.seed(1)
#' n <- 50
#' m <- 30
#' argvals <- seq(0, 1, length.out = m)
#' X <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(X, argvals = argvals)
#'
#' # Build Phase I chart
#' chart <- spm.phase1(fd, ncomp = 3)
#' chart
#' }
#'
#' @export
spm.phase1 <- function(fdataobj, ncomp = 5, alpha = 0.05,
                        tuning.fraction = 0.5, seed = 42) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  n <- nrow(fdataobj$data)
  m <- ncol(fdataobj$data)

  if (n < 4) {
    stop("At least 4 observations are required for SPM Phase I")
  }

  argvals <- as.numeric(fdataobj$argvals)

  result <- .Call("wrap__spm_phase1_rust",
                  fdataobj$data, argvals,
                  as.integer(ncomp), as.numeric(alpha),
                  as.numeric(tuning.fraction), as.integer(seed))

  if (is.null(result)) {
    stop("spm.phase1 failed: check data dimensions and parameters")
  }

  structure(
    list(
      eigenvalues = result$eigenvalues,
      t2.phase1 = result$t2_phase1,
      spe.phase1 = result$spe_phase1,
      t2.ucl = result$t2_ucl,
      spe.ucl = result$spe_ucl,
      ncomp = result$ncomp,
      alpha = alpha,
      t2.description = result$t2_description,
      spe.description = result$spe_description,
      fdataobj = fdataobj,
      call = match.call(),
      .rust = list(
        rotation = result$rotation,
        scores = result$scores,
        mean = result$mean,
        singular_values = result$singular_values,
        centered = result$centered,
        eigenvalues = result$eigenvalues,
        t2_ucl = result$t2_ucl,
        t2_alpha = result$t2_alpha,
        t2_description = result$t2_description,
        spe_ucl = result$spe_ucl,
        spe_alpha = result$spe_alpha,
        spe_description = result$spe_description,
        ncomp = result$ncomp,
        config_alpha = result$config_alpha,
        config_tuning_fraction = result$config_tuning_fraction,
        config_seed = result$config_seed
      )
    ),
    class = "spm.chart"
  )
}

# =============================================================================
# Phase II: Monitor new data
# =============================================================================

#' Monitor New Functional Data (Phase II)
#'
#' Projects new functional observations through a trained SPM chart and
#' checks whether T-squared and SPE statistics exceed control limits.
#'
#' @param chart An object of class \code{spm.chart} from \code{\link{spm.phase1}}.
#' @param newdata An object of class \code{fdata} with new observations to monitor.
#'
#' @return An object of class \code{spm.monitor} with components:
#' \describe{
#'   \item{t2}{T-squared values for new observations}
#'   \item{spe}{SPE values for new observations}
#'   \item{t2.alarm}{Logical vector: TRUE where T-squared exceeds UCL}
#'   \item{spe.alarm}{Logical vector: TRUE where SPE exceeds UCL}
#'   \item{scores}{FPC score matrix for new observations}
#'   \item{t2.ucl}{T-squared control limit}
#'   \item{spe.ucl}{SPE control limit}
#' }
#'
#' @seealso \code{\link{spm.phase1}} for building the chart
#'
#' @examples
#' \donttest{
#' # Build chart from in-control data
#' set.seed(1)
#' n <- 50; m <- 30
#' argvals <- seq(0, 1, length.out = m)
#' X <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(X, argvals = argvals)
#' chart <- spm.phase1(fd, ncomp = 3)
#'
#' # Monitor new data (with a shift)
#' X_new <- matrix(rnorm(10 * m) + 2, 10, m)
#' fd_new <- fdata(X_new, argvals = argvals)
#' mon <- spm.monitor(chart, fd_new)
#' mon
#' }
#'
#' @export
spm.monitor <- function(chart, newdata) {
  if (!inherits(chart, "spm.chart")) {
    stop("chart must be of class 'spm.chart'")
  }
  if (!inherits(newdata, "fdata")) {
    stop("newdata must be of class 'fdata'")
  }

  r <- chart$.rust
  argvals <- as.numeric(chart$fdataobj$argvals)

  result <- .Call("wrap__spm_monitor_rust",
                  newdata$data, argvals,
                  r$rotation, r$mean, r$singular_values,
                  r$centered,
                  r$eigenvalues,
                  r$t2_ucl, r$t2_alpha, r$t2_description,
                  r$spe_ucl, r$spe_alpha, r$spe_description,
                  as.integer(r$ncomp), r$config_alpha,
                  r$config_tuning_fraction, as.integer(r$config_seed))

  if (is.null(result)) {
    stop("spm.monitor failed: check new data dimensions")
  }

  structure(
    list(
      t2 = result$t2,
      spe = result$spe,
      t2.alarm = as.logical(result$t2_alarm),
      spe.alarm = as.logical(result$spe_alarm),
      scores = result$scores,
      t2.ucl = chart$t2.ucl,
      spe.ucl = chart$spe.ucl,
      call = match.call()
    ),
    class = "spm.monitor"
  )
}

# =============================================================================
# MFPCA: Multivariate Functional PCA
# =============================================================================

#' Multivariate Functional Principal Component Analysis
#'
#' Extends univariate FPCA to handle multiple functional variables
#' observed on potentially different grids. Variables are optionally
#' weighted by their inverse standard deviation before joint SVD.
#'
#' @param fdataobj.list A list of \code{fdata} objects, one per functional variable.
#'   All must have the same number of observations (rows) but may differ
#'   in grid size.
#' @param ncomp Number of principal components to extract (default 5).
#' @param weighted Logical; whether to weight each variable by 1/std_dev
#'   before SVD (default TRUE).
#'
#' @return An object of class \code{mfpca} with components:
#' \describe{
#'   \item{scores}{Score matrix (n x ncomp)}
#'   \item{eigenfunctions}{List of eigenfunction matrices, one per variable}
#'   \item{eigenvalues}{Eigenvalues (length ncomp)}
#'   \item{means}{List of mean functions, one per variable}
#'   \item{scales}{Per-variable standard deviations}
#'   \item{grid.sizes}{Grid sizes per variable}
#' }
#'
#' @examples
#' \donttest{
#' # Two functional variables
#' n <- 40
#' fd1 <- fdata(matrix(rnorm(n * 20), n, 20), argvals = seq(0, 1, length.out = 20))
#' fd2 <- fdata(matrix(rnorm(n * 15), n, 15), argvals = seq(0, 1, length.out = 15))
#' result <- mfpca(list(fd1, fd2), ncomp = 3)
#' dim(result$scores)
#' }
#'
#' @export
mfpca <- function(fdataobj.list, ncomp = 5, weighted = TRUE) {
  if (!is.list(fdataobj.list)) {
    stop("fdataobj.list must be a list of fdata objects")
  }
  if (length(fdataobj.list) == 0) {
    stop("fdataobj.list must contain at least one fdata object")
  }

  for (i in seq_along(fdataobj.list)) {
    if (!inherits(fdataobj.list[[i]], "fdata")) {
      stop("Element ", i, " of fdataobj.list is not an fdata object")
    }
  }

  # Extract data matrices
  data_list <- lapply(fdataobj.list, function(fd) fd$data)

  result <- .Call("wrap__mfpca_rust", data_list,
                  as.integer(ncomp), as.logical(weighted))

  if (is.null(result)) {
    stop("mfpca failed: check data dimensions")
  }

  structure(
    list(
      scores = result$scores,
      eigenfunctions = result$eigenfunctions,
      eigenvalues = result$eigenvalues,
      means = result$means,
      scales = result$scales,
      grid.sizes = result$grid_sizes,
      fdataobj.list = fdataobj.list,
      call = match.call(),
      .rust = list()
    ),
    class = "mfpca"
  )
}

# =============================================================================
# EWMA monitoring
# =============================================================================

#' EWMA-Based SPM Monitoring
#'
#' Applies Exponentially Weighted Moving Average (EWMA) smoothing to FPC
#' scores before computing monitoring statistics. This increases sensitivity
#' to small persistent shifts in the process.
#'
#' @param chart An object of class \code{spm.chart} from \code{\link{spm.phase1}}.
#' @param newdata An object of class \code{fdata} with sequential observations
#'   (rows in time order).
#' @param lambda EWMA smoothing parameter in (0, 1] (default 0.2).
#'   Smaller values give more smoothing; lambda = 1 gives raw scores.
#' @param ncomp Number of components for EWMA (default: same as chart).
#' @param alpha Significance level for EWMA control limit (default 0.05).
#'
#' @return A list with components:
#' \describe{
#'   \item{smoothed.scores}{EWMA-smoothed score matrix}
#'   \item{t2}{T-squared values on smoothed scores}
#'   \item{spe}{SPE values (unsmoothed)}
#'   \item{t2.limit}{T-squared control limit for EWMA}
#'   \item{spe.limit}{SPE control limit}
#'   \item{t2.alarm}{Logical: TRUE where T-squared exceeds limit}
#'   \item{spe.alarm}{Logical: TRUE where SPE exceeds limit}
#' }
#'
#' @seealso \code{\link{spm.phase1}} for building the chart,
#'   \code{\link{spm.monitor}} for standard (non-EWMA) monitoring
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 50; m <- 30
#' argvals <- seq(0, 1, length.out = m)
#' X <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(X, argvals = argvals)
#' chart <- spm.phase1(fd, ncomp = 3)
#'
#' # Sequential monitoring with small shift
#' X_seq <- matrix(rnorm(20 * m) + 0.5, 20, m)
#' fd_seq <- fdata(X_seq, argvals = argvals)
#' ewma_result <- spm.ewma(chart, fd_seq, lambda = 0.2)
#' ewma_result$t2.alarm
#' }
#'
#' @export
spm.ewma <- function(chart, newdata, lambda = 0.2, ncomp = NULL, alpha = 0.05) {
  if (!inherits(chart, "spm.chart")) {
    stop("chart must be of class 'spm.chart'")
  }
  if (!inherits(newdata, "fdata")) {
    stop("newdata must be of class 'fdata'")
  }

  if (is.null(ncomp)) {
    ncomp <- chart$ncomp
  }

  r <- chart$.rust
  argvals <- as.numeric(chart$fdataobj$argvals)

  result <- .Call("wrap__spm_ewma_rust",
                  newdata$data, argvals,
                  r$rotation, r$mean, r$singular_values,
                  r$centered,
                  r$eigenvalues,
                  r$t2_ucl, r$t2_alpha, r$t2_description,
                  r$spe_ucl, r$spe_alpha, r$spe_description,
                  as.integer(r$ncomp), r$config_alpha,
                  r$config_tuning_fraction, as.integer(r$config_seed),
                  as.numeric(lambda), as.integer(ncomp), as.numeric(alpha))

  if (is.null(result)) {
    stop("spm.ewma failed: check data dimensions and parameters")
  }

  list(
    smoothed.scores = result$smoothed_scores,
    t2 = result$t2,
    spe = result$spe,
    t2.limit = result$t2_limit,
    spe.limit = result$spe_limit,
    t2.alarm = as.logical(result$t2_alarm),
    spe.alarm = as.logical(result$spe_alarm)
  )
}

# =============================================================================
# Contribution diagnostics
# =============================================================================

#' SPM Contribution Diagnostics
#'
#' When an alarm is triggered, identifies which functional variables or
#' components contribute most to the elevated T-squared or SPE statistic.
#' Useful for root-cause analysis.
#'
#' @param scores Score matrix (n x ncomp).
#' @param eigenvalues Eigenvalues vector (length ncomp). Required for T-squared.
#' @param grid.sizes Integer vector of grid sizes per variable. For univariate
#'   SPM, use a single value equal to ncomp.
#' @param type Character; either \code{"t2"} for T-squared contributions or
#'   \code{"spe"} for SPE contributions. Default \code{"t2"}.
#' @param standardized.vars For SPE contributions: list of standardized data matrices.
#' @param reconstructed.vars For SPE contributions: list of reconstructed data matrices.
#' @param argvals.list For SPE contributions: list of argvals vectors.
#'
#' @return A matrix (n x p) of per-variable contributions.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' scores <- matrix(rnorm(10 * 3), 10, 3)
#' eigenvalues <- c(5, 2, 1)
#' # Per-component contributions (univariate, 1 variable with 3 grid points)
#' contrib <- spm.contributions(scores, eigenvalues, grid.sizes = 3L)
#' dim(contrib)
#' }
#'
#' @export
spm.contributions <- function(scores, eigenvalues = NULL, grid.sizes = NULL,
                               type = c("t2", "spe"),
                               standardized.vars = NULL,
                               reconstructed.vars = NULL,
                               argvals.list = NULL) {
  type <- match.arg(type)

  if (type == "t2") {
    if (is.null(eigenvalues)) {
      stop("eigenvalues required for T-squared contributions")
    }
    if (is.null(grid.sizes)) {
      # Default: one variable with ncomp grid points
      grid.sizes <- ncol(scores)
    }

    result <- .Call("wrap__spm_t2_contrib_rust",
                    scores, as.numeric(eigenvalues), as.integer(grid.sizes))

    if (is.null(result)) {
      stop("T-squared contribution computation failed")
    }
    return(result)
  } else {
    # SPE contributions
    if (is.null(standardized.vars) || is.null(reconstructed.vars) ||
        is.null(argvals.list)) {
      stop("standardized.vars, reconstructed.vars, and argvals.list required for SPE contributions")
    }

    result <- .Call("wrap__spm_spe_contrib_rust",
                    standardized.vars, reconstructed.vars, argvals.list)

    if (is.null(result)) {
      stop("SPE contribution computation failed")
    }
    return(result)
  }
}

# =============================================================================
# FRCC: Functional Regression Control Chart
# =============================================================================

#' Build Functional Regression Control Chart (Phase I)
#'
#' Monitors a functional response after adjusting for known scalar covariates
#' using function-on-scalar regression (FOSR). The residuals are then monitored
#' via FPCA-based T-squared and SPE statistics.
#'
#' @param fdataobj An object of class \code{fdata} (functional response).
#' @param predictors A matrix of scalar predictors (n x p).
#' @param ncomp Number of principal components for residual FPCA (default 5).
#' @param fosr.lambda FOSR smoothing parameter (default 1e-4).
#' @param alpha Significance level (default 0.05).
#' @param tuning.fraction Fraction of data for tuning (default 0.5).
#' @param seed Random seed (default 42).
#'
#' @return An object of class \code{frcc.chart} with components:
#' \describe{
#'   \item{eigenvalues}{Eigenvalues from residual FPCA}
#'   \item{t2.ucl}{T-squared control limit}
#'   \item{spe.ucl}{SPE control limit}
#'   \item{ncomp}{Number of components used}
#'   \item{fdataobj}{Original fdata object}
#'   \item{predictors}{Original predictor matrix}
#'   \item{.rust}{Internal fields for Phase II monitoring}
#' }
#'
#' @seealso \code{\link{frcc.monitor}} for Phase II monitoring,
#'   \code{\link{spm.phase1}} for monitoring without covariates
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 60; m <- 20
#' argvals <- seq(0, 1, length.out = m)
#' X_pred <- cbind(rnorm(n), rnorm(n))
#' Y <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(Y, argvals = argvals)
#'
#' chart <- frcc.phase1(fd, X_pred, ncomp = 3)
#' chart
#' }
#'
#' @export
frcc.phase1 <- function(fdataobj, predictors, ncomp = 5, fosr.lambda = 1e-4,
                         alpha = 0.05, tuning.fraction = 0.5, seed = 42) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  predictors <- as.matrix(predictors)
  n <- nrow(fdataobj$data)

  if (nrow(predictors) != n) {
    stop("Number of rows in predictors must equal number of curves")
  }

  argvals <- as.numeric(fdataobj$argvals)

  result <- .Call("wrap__frcc_phase1_rust",
                  fdataobj$data, predictors, argvals,
                  as.integer(ncomp), as.numeric(fosr.lambda),
                  as.numeric(alpha), as.numeric(tuning.fraction),
                  as.integer(seed))

  if (is.null(result)) {
    stop("frcc.phase1 failed: check data dimensions and parameters")
  }

  structure(
    list(
      eigenvalues = result$eigenvalues,
      t2.ucl = result$t2_ucl,
      spe.ucl = result$spe_ucl,
      ncomp = result$ncomp,
      alpha = alpha,
      t2.description = result$t2_description,
      spe.description = result$spe_description,
      fdataobj = fdataobj,
      predictors = predictors,
      call = match.call(),
      .rust = list(
        fosr_intercept = result$fosr_intercept,
        fosr_beta = result$fosr_beta,
        fosr_lambda = result$fosr_lambda,
        resid_rotation = result$resid_rotation,
        resid_mean = result$resid_mean,
        resid_singular_values = result$resid_singular_values,
        resid_centered = result$resid_centered,
        eigenvalues = result$eigenvalues,
        t2_ucl = result$t2_ucl,
        t2_alpha = result$t2_alpha,
        t2_description = result$t2_description,
        spe_ucl = result$spe_ucl,
        spe_alpha = result$spe_alpha,
        spe_description = result$spe_description,
        ncomp = result$ncomp,
        config_fosr_lambda = result$config_fosr_lambda,
        config_alpha = result$config_alpha,
        config_tuning_fraction = result$config_tuning_fraction,
        config_seed = result$config_seed
      )
    ),
    class = "frcc.chart"
  )
}

#' Monitor New Data Against FRCC Chart (Phase II)
#'
#' Predicts the functional response from the FOSR model, computes
#' residuals, and monitors them against the established FRCC chart.
#'
#' @param chart An object of class \code{frcc.chart} from \code{\link{frcc.phase1}}.
#' @param newdata An object of class \code{fdata} with new functional response.
#' @param new.predictors A matrix of new scalar predictors (n_new x p).
#'
#' @return An object of class \code{spm.monitor} with components:
#' \describe{
#'   \item{t2}{T-squared values}
#'   \item{spe}{SPE values}
#'   \item{t2.alarm}{Logical: TRUE where T-squared exceeds UCL}
#'   \item{spe.alarm}{Logical: TRUE where SPE exceeds UCL}
#'   \item{residual.scores}{Residual FPC scores}
#'   \item{t2.ucl}{T-squared control limit}
#'   \item{spe.ucl}{SPE control limit}
#' }
#'
#' @seealso \code{\link{frcc.phase1}} for building the chart
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 60; m <- 20
#' argvals <- seq(0, 1, length.out = m)
#' X_pred <- cbind(rnorm(n), rnorm(n))
#' Y <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(Y, argvals = argvals)
#' chart <- frcc.phase1(fd, X_pred, ncomp = 3)
#'
#' # Monitor new data
#' X_new <- cbind(rnorm(10), rnorm(10))
#' Y_new <- matrix(rnorm(10 * m) + 2, 10, m)
#' fd_new <- fdata(Y_new, argvals = argvals)
#' mon <- frcc.monitor(chart, fd_new, X_new)
#' mon
#' }
#'
#' @export
frcc.monitor <- function(chart, newdata, new.predictors) {
  if (!inherits(chart, "frcc.chart")) {
    stop("chart must be of class 'frcc.chart'")
  }
  if (!inherits(newdata, "fdata")) {
    stop("newdata must be of class 'fdata'")
  }

  new.predictors <- as.matrix(new.predictors)
  r <- chart$.rust
  argvals <- as.numeric(chart$fdataobj$argvals)

  result <- .Call("wrap__frcc_monitor_rust",
                  newdata$data, new.predictors, argvals,
                  r$fosr_intercept, r$fosr_beta,
                  r$fosr_lambda,
                  r$resid_rotation, r$resid_mean,
                  r$resid_singular_values, r$resid_centered,
                  r$eigenvalues,
                  r$t2_ucl, r$t2_alpha, r$t2_description,
                  r$spe_ucl, r$spe_alpha, r$spe_description,
                  as.integer(r$ncomp), r$config_fosr_lambda,
                  r$config_alpha, r$config_tuning_fraction,
                  as.integer(r$config_seed))

  if (is.null(result)) {
    stop("frcc.monitor failed: check data dimensions")
  }

  structure(
    list(
      t2 = result$t2,
      spe = result$spe,
      t2.alarm = as.logical(result$t2_alarm),
      spe.alarm = as.logical(result$spe_alarm),
      residual.scores = result$residual_scores,
      t2.ucl = chart$t2.ucl,
      spe.ucl = chart$spe.ucl,
      call = match.call()
    ),
    class = "spm.monitor"
  )
}

# =============================================================================
# Plot methods
# =============================================================================

#' Plot an SPM Phase I control chart
#'
#' Displays T-squared and SPE statistics for Phase I observations with
#' upper control limits.
#'
#' @param x An object of class \code{spm.chart}.
#' @param ... Additional arguments (ignored).
#' @export
plot.spm.chart <- function(x, ...) {
  n <- length(x$t2.phase1)
  df <- data.frame(
    obs = rep(seq_len(n), 2),
    value = c(x$t2.phase1, x$spe.phase1),
    panel = factor(
      rep(c("Hotelling T\u00b2", "Squared Prediction Error"), each = n),
      levels = c("Hotelling T\u00b2", "Squared Prediction Error")
    ),
    ucl = rep(c(x$t2.ucl, x$spe.ucl), each = n)
  )
  df$alarm <- df$value > df$ucl

  ucl_df <- data.frame(
    panel = factor(
      c("Hotelling T\u00b2", "Squared Prediction Error"),
      levels = c("Hotelling T\u00b2", "Squared Prediction Error")
    ),
    ucl = c(x$t2.ucl, x$spe.ucl),
    label = paste("UCL =", format(c(x$t2.ucl, x$spe.ucl), digits = 4))
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$obs, y = .data$value)) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = .data$obs, y = 0, yend = .data$value,
                   color = .data$alarm),
      linewidth = 0.8
    ) +
    ggplot2::geom_hline(
      data = ucl_df,
      ggplot2::aes(yintercept = .data$ucl),
      linetype = "dashed", color = "#D55E00", linewidth = 0.7
    ) +
    ggplot2::geom_text(
      data = ucl_df,
      ggplot2::aes(x = -Inf, y = .data$ucl, label = .data$label),
      hjust = -0.05, vjust = -0.5, size = 3, color = "#D55E00"
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "#4A90D9", "TRUE" = "#D55E00"),
      guide = "none"
    ) +
    ggplot2::facet_wrap(~ .data$panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = "Phase I Control Chart",
      x = "Observation", y = "Statistic"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))

  p
}

#' Plot SPM Phase II monitoring results
#'
#' Displays T-squared and SPE monitoring statistics with upper control limits.
#' Alarm observations are highlighted in red.
#'
#' @param x An object of class \code{spm.monitor}.
#' @param ... Additional arguments (ignored).
#' @export
plot.spm.monitor <- function(x, ...) {
  n <- length(x$t2)
  df <- data.frame(
    obs = rep(seq_len(n), 2),
    value = c(x$t2, x$spe),
    panel = factor(
      rep(c("Hotelling T\u00b2", "Squared Prediction Error"), each = n),
      levels = c("Hotelling T\u00b2", "Squared Prediction Error")
    ),
    ucl = rep(c(x$t2.ucl, x$spe.ucl), each = n),
    alarm = c(x$t2.alarm, x$spe.alarm)
  )

  ucl_df <- data.frame(
    panel = factor(
      c("Hotelling T\u00b2", "Squared Prediction Error"),
      levels = c("Hotelling T\u00b2", "Squared Prediction Error")
    ),
    ucl = c(x$t2.ucl, x$spe.ucl),
    label = paste("UCL =", format(c(x$t2.ucl, x$spe.ucl), digits = 4))
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$obs, y = .data$value)) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = .data$obs, y = 0, yend = .data$value,
                   color = .data$alarm),
      linewidth = 0.8
    ) +
    ggplot2::geom_hline(
      data = ucl_df,
      ggplot2::aes(yintercept = .data$ucl),
      linetype = "dashed", color = "#D55E00", linewidth = 0.7
    ) +
    ggplot2::geom_text(
      data = ucl_df,
      ggplot2::aes(x = -Inf, y = .data$ucl, label = .data$label),
      hjust = -0.05, vjust = -0.5, size = 3, color = "#D55E00"
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "#4A90D9", "TRUE" = "#D55E00"),
      guide = "none"
    ) +
    ggplot2::facet_wrap(~ .data$panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = "Phase II Monitoring",
      x = "Observation", y = "Statistic"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))

  p
}

# =============================================================================
# Print methods
# =============================================================================

#' @export
print.spm.chart <- function(x, ...) {
  cat("SPM Control Chart (Phase I)\n")
  cat("  Components:", x$ncomp, "\n")
  cat("  Alpha:", x$alpha, "\n")
  cat("  T2 UCL:", format(x$t2.ucl, digits = 4), "\n")
  cat("  SPE UCL:", format(x$spe.ucl, digits = 4), "\n")
  cat("  Observations:", nrow(x$fdataobj$data), "\n")
  cat("  Grid points:", ncol(x$fdataobj$data), "\n")
  cat("  Eigenvalues:", paste(format(x$eigenvalues, digits = 3), collapse = ", "), "\n")
  invisible(x)
}

#' @export
print.spm.monitor <- function(x, ...) {
  n <- length(x$t2)
  n_t2 <- sum(x$t2.alarm)
  n_spe <- sum(x$spe.alarm)

  cat("SPM Monitoring Result (Phase II)\n")
  cat("  Observations:", n, "\n")
  cat("  T2 alarms:", n_t2, "of", n,
      paste0("(", format(100 * n_t2 / n, digits = 3), "%)"), "\n")
  cat("  SPE alarms:", n_spe, "of", n,
      paste0("(", format(100 * n_spe / n, digits = 3), "%)"), "\n")
  cat("  T2 UCL:", format(x$t2.ucl, digits = 4), "\n")
  cat("  SPE UCL:", format(x$spe.ucl, digits = 4), "\n")
  invisible(x)
}

#' @export
print.mfpca <- function(x, ...) {
  cat("Multivariate Functional PCA\n")
  cat("  Variables:", length(x$grid.sizes), "\n")
  cat("  Grid sizes:", paste(x$grid.sizes, collapse = ", "), "\n")
  cat("  Components:", length(x$eigenvalues), "\n")
  pve <- x$eigenvalues / sum(x$eigenvalues) * 100
  cat("  Variance explained:", paste(format(pve, digits = 3), collapse = ", "), "%\n")
  invisible(x)
}

#' @export
print.frcc.chart <- function(x, ...) {
  cat("Functional Regression Control Chart (Phase I)\n")
  cat("  Components:", x$ncomp, "\n")
  cat("  Alpha:", x$alpha, "\n")
  cat("  T2 UCL:", format(x$t2.ucl, digits = 4), "\n")
  cat("  SPE UCL:", format(x$spe.ucl, digits = 4), "\n")
  cat("  Observations:", nrow(x$fdataobj$data), "\n")
  cat("  Predictors:", ncol(x$predictors), "\n")
  invisible(x)
}

# =============================================================================
# Phase 1 Extensions: ncomp selection, PC contributions, rules, robust limits
# =============================================================================

#' Select Number of Principal Components
#'
#' Selects the number of principal components to retain based on eigenvalues
#' using one of several criteria: cumulative variance threshold, Kaiser rule,
#' elbow method, or a fixed number.
#'
#' @param eigenvalues Numeric vector of eigenvalues from FPCA.
#' @param method Selection method. One of \code{"variance90"} (cumulative
#'   variance exceeds \code{threshold}), \code{"kaiser"} (eigenvalues > mean),
#'   \code{"elbow"} (scree plot elbow), or \code{"fixed"} (use \code{threshold}
#'   as the number of components).
#' @param threshold For \code{"variance90"}, the cumulative proportion of
#'   variance to exceed (default 0.9). For \code{"fixed"}, the exact number
#'   of components.
#'
#' @return An integer: the selected number of components.
#'
#' @examples
#' \donttest{
#' eigenvalues <- c(5.0, 2.0, 1.0, 0.5, 0.2, 0.1)
#' spm.ncomp.select(eigenvalues, method = "variance90", threshold = 0.9)
#' spm.ncomp.select(eigenvalues, method = "kaiser")
#' }
#'
#' @export
spm.ncomp.select <- function(eigenvalues,
                              method = c("variance90", "kaiser", "elbow", "fixed"),
                              threshold = 0.9) {
  method <- match.arg(method)
  eigenvalues <- as.numeric(eigenvalues)

  if (length(eigenvalues) == 0) {
    stop("eigenvalues must be a non-empty numeric vector")
  }

  result <- .Call("wrap__spm_select_ncomp_rust",
                  eigenvalues, method, as.numeric(threshold))

  if (is.null(result)) {
    stop("spm.ncomp.select failed: check eigenvalues and method")
  }

  as.integer(result)
}

#' Per-PC Contributions to T-Squared Statistic
#'
#' Computes the contribution of each principal component to the Hotelling
#' T-squared statistic for each observation. Useful for diagnosing which
#' components drive an alarm.
#'
#' @param scores Score matrix (n x ncomp), e.g., from \code{spm.monitor()$scores}.
#' @param eigenvalues Numeric vector of eigenvalues (length ncomp).
#'
#' @return A matrix (n x ncomp) of per-component T-squared contributions.
#'   Row sums equal the T-squared statistics.
#'
#' @seealso \code{\link{spm.contributions}} for per-variable contributions
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' scores <- matrix(rnorm(10 * 3), 10, 3)
#' eigenvalues <- c(5, 2, 1)
#' contrib <- spm.pc.contributions(scores, eigenvalues)
#' dim(contrib)  # 10 x 3
#' }
#'
#' @export
spm.pc.contributions <- function(scores, eigenvalues) {
  scores <- as.matrix(scores)
  eigenvalues <- as.numeric(eigenvalues)

  if (ncol(scores) != length(eigenvalues)) {
    stop("Number of columns in scores must equal length of eigenvalues")
  }

  result <- .Call("wrap__spm_t2_pc_contrib_rust",
                  scores, eigenvalues)

  if (is.null(result)) {
    stop("spm.pc.contributions failed: check scores and eigenvalues")
  }

  result
}

#' Evaluate Control Chart Rules
#'
#' Applies Western Electric or Nelson run rules to a sequence of control
#' chart values. Identifies patterns that indicate non-random behavior
#' (e.g., runs, trends, points beyond control limits).
#'
#' @param values Numeric vector of chart values (in time order).
#' @param center Center line value (e.g., process mean).
#' @param sigma Standard deviation (for zone boundaries at 1, 2, 3 sigma).
#' @param rule.set Which rule set to apply: \code{"western.electric"} (4 rules)
#'   or \code{"nelson"} (8 rules).
#'
#' @return An object of class \code{spm.rules} with components:
#' \describe{
#'   \item{rules}{Character vector of violated rule names}
#'   \item{indices}{List of integer vectors; each gives the observation
#'     indices (1-indexed) involved in that violation}
#'   \item{values}{The input values}
#'   \item{center}{The center line}
#'   \item{sigma}{The standard deviation}
#'   \item{rule.set}{The rule set used}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' values <- c(rnorm(20), rnorm(5, mean = 3))
#' result <- spm.rules(values, center = 0, sigma = 1,
#'                     rule.set = "western.electric")
#' result
#' }
#'
#' @export
spm.rules <- function(values, center, sigma,
                       rule.set = c("western.electric", "nelson")) {
  rule.set <- match.arg(rule.set)
  values <- as.numeric(values)
  center <- as.numeric(center)
  sigma <- as.numeric(sigma)

  if (length(values) == 0) {
    stop("values must be a non-empty numeric vector")
  }
  if (sigma <= 0) {
    stop("sigma must be positive")
  }

  # Map R-style dots to Rust-style underscores
  rule_set_rust <- gsub("\\.", "_", rule.set)

  result <- .Call("wrap__spm_evaluate_rules_rust",
                  values, center, sigma, rule_set_rust)

  if (is.null(result)) {
    stop("spm.rules failed: check input values")
  }

  structure(
    list(
      rules = result$rules,
      indices = result$indices,
      values = values,
      center = center,
      sigma = sigma,
      rule.set = rule.set,
      call = match.call()
    ),
    class = "spm.rules"
  )
}

#' @export
print.spm.rules <- function(x, ...) {
  cat("SPM Control Chart Rule Evaluation\n")
  cat("  Rule set:", x$rule.set, "\n")
  cat("  Observations:", length(x$values), "\n")
  cat("  Center:", format(x$center, digits = 4), "\n")
  cat("  Sigma:", format(x$sigma, digits = 4), "\n")
  n_violations <- length(x$rules)
  cat("  Violations found:", n_violations, "\n")
  if (n_violations > 0) {
    for (i in seq_len(n_violations)) {
      idx <- x$indices[[i]]
      idx_str <- if (length(idx) > 10) {
        paste(c(paste(idx[1:10], collapse = ", "), "..."), collapse = "")
      } else {
        paste(idx, collapse = ", ")
      }
      cat("    ", x$rules[i], ": indices [", idx_str, "]\n", sep = "")
    }
  }
  invisible(x)
}

#' Robust Control Limits via Alternative Methods
#'
#' Computes T-squared or SPE control limits using methods beyond the standard
#' parametric approach: empirical quantile, bootstrap, or kernel density
#' estimation. Useful when the parametric chi-squared assumption is questionable.
#'
#' @param values Numeric vector of Phase I statistic values (T-squared or SPE).
#' @param ncomp Number of components (required for T-squared; ignored for SPE).
#' @param alpha Significance level (default 0.05).
#' @param type Which statistic: \code{"t2"} or \code{"spe"}.
#' @param method Estimation method: \code{"parametric"}, \code{"empirical"},
#'   \code{"bootstrap"}, or \code{"kde"}.
#' @param n.bootstrap Number of bootstrap resamples (default 1000; only used
#'   when \code{method = "bootstrap"}).
#' @param seed Random seed for bootstrap (default 42).
#'
#' @return A list with components:
#' \describe{
#'   \item{ucl}{Upper control limit}
#'   \item{alpha}{Significance level used}
#'   \item{description}{Character string describing the method}
#'   \item{method}{The method used}
#' }
#'
#' @seealso \code{\link{spm.phase1}} for standard control limits
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' t2_values <- rchisq(100, df = 3)
#' spm.limit.robust(t2_values, ncomp = 3, type = "t2", method = "empirical")
#' spm.limit.robust(t2_values, ncomp = 3, type = "t2", method = "bootstrap")
#' }
#'
#' @export
spm.limit.robust <- function(values, ncomp = NULL, alpha = 0.05,
                              type = c("t2", "spe"),
                              method = c("parametric", "empirical",
                                         "bootstrap", "kde"),
                              n.bootstrap = 1000, seed = 42) {
  type <- match.arg(type)
  method <- match.arg(method)
  values <- as.numeric(values)

  if (length(values) == 0) {
    stop("values must be a non-empty numeric vector")
  }

  if (type == "t2") {
    if (is.null(ncomp)) {
      stop("ncomp is required for T-squared limits")
    }
    result <- .Call("wrap__spm_t2_limit_robust_rust",
                    values, as.integer(ncomp), as.numeric(alpha),
                    method, as.integer(n.bootstrap), as.integer(seed))
  } else {
    result <- .Call("wrap__spm_spe_limit_robust_rust",
                    values, as.numeric(alpha),
                    method, as.integer(n.bootstrap), as.integer(seed))
  }

  if (is.null(result)) {
    stop("spm.limit.robust failed: check input values and parameters")
  }

  list(
    ucl = result$ucl,
    alpha = result$alpha,
    description = result$description,
    method = method
  )
}

# =============================================================================
# Phase 2: ARL estimation and CUSUM monitoring
# =============================================================================

#' Estimate Average Run Length (ARL)
#'
#' Estimates the Average Run Length (ARL) for a T-squared or SPE control chart
#' via Monte Carlo simulation. In-control ARL (\code{shift.magnitude = NULL})
#' measures the expected number of observations before a false alarm.
#' Out-of-control ARL measures detection speed under a mean shift.
#'
#' @param chart An object of class \code{spm.chart} from \code{\link{spm.phase1}}.
#' @param type Which statistic: \code{"t2"} or \code{"spe"}.
#' @param shift.magnitude For out-of-control ARL: numeric vector of per-component
#'   shifts (length = ncomp). If \code{NULL}, computes in-control ARL.
#' @param n.sim Number of Monte Carlo simulations (default 10000).
#' @param max.rl Maximum run length per simulation (default 5000).
#' @param seed Random seed (default 42).
#'
#' @return An object of class \code{spm.arl} with components:
#' \describe{
#'   \item{arl}{Estimated ARL}
#'   \item{std.dev}{Standard deviation of run lengths}
#'   \item{median.rl}{Median run length}
#'   \item{run.lengths}{Integer vector of all simulated run lengths}
#'   \item{type}{Which statistic was used}
#'   \item{in.control}{Logical: TRUE if in-control ARL}
#' }
#'
#' @seealso \code{\link{spm.arl.ewma}} for EWMA-based ARL
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 50; m <- 30
#' argvals <- seq(0, 1, length.out = m)
#' X <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(X, argvals = argvals)
#' chart <- spm.phase1(fd, ncomp = 3)
#'
#' # In-control ARL
#' arl0 <- spm.arl(chart, type = "t2", n.sim = 1000)
#' arl0
#' }
#'
#' @export
spm.arl <- function(chart, type = c("t2", "spe"),
                     shift.magnitude = NULL,
                     n.sim = 10000, max.rl = 5000, seed = 42) {
  if (!inherits(chart, "spm.chart")) {
    stop("chart must be of class 'spm.chart'")
  }
  type <- match.arg(type)

  if (type == "t2") {
    if (is.null(shift.magnitude)) {
      result <- .Call("wrap__spm_arl0_t2_rust",
                      as.numeric(chart$eigenvalues),
                      as.numeric(chart$t2.ucl),
                      as.integer(n.sim), as.integer(max.rl),
                      as.integer(seed))
    } else {
      shift.magnitude <- as.numeric(shift.magnitude)
      if (length(shift.magnitude) != chart$ncomp) {
        stop("shift.magnitude must have length equal to ncomp (", chart$ncomp, ")")
      }
      result <- .Call("wrap__spm_arl1_t2_rust",
                      as.numeric(chart$eigenvalues),
                      as.numeric(chart$t2.ucl),
                      shift.magnitude,
                      as.integer(n.sim), as.integer(max.rl),
                      as.integer(seed))
    }
  } else {
    # SPE ARL requires df and scale parameters
    # Estimate from Phase I SPE values using method of moments
    spe_vals <- chart$spe.phase1
    spe_mean <- mean(spe_vals)
    spe_var <- stats::var(spe_vals)
    spe_scale <- spe_var / (2 * spe_mean)
    spe_df <- 2 * spe_mean^2 / spe_var

    result <- .Call("wrap__spm_arl0_spe_rust",
                    as.numeric(spe_df), as.numeric(spe_scale),
                    as.numeric(chart$spe.ucl),
                    as.integer(n.sim), as.integer(max.rl),
                    as.integer(seed))
  }

  if (is.null(result)) {
    stop("spm.arl failed: check chart and parameters")
  }

  structure(
    list(
      arl = result$arl,
      std.dev = result$std_dev,
      median.rl = result$median_rl,
      run.lengths = result$run_lengths,
      type = type,
      in.control = is.null(shift.magnitude),
      call = match.call()
    ),
    class = "spm.arl"
  )
}

#' @export
print.spm.arl <- function(x, ...) {
  label <- if (x$in.control) "In-Control" else "Out-of-Control"
  cat("SPM Average Run Length (", label, ")\n", sep = "")
  cat("  Statistic:", toupper(x$type), "\n")
  cat("  ARL:", format(x$arl, digits = 4), "\n")
  cat("  Std. Dev.:", format(x$std.dev, digits = 4), "\n")
  cat("  Median RL:", x$median.rl, "\n")
  cat("  Simulations:", length(x$run.lengths), "\n")
  invisible(x)
}

#' Estimate ARL for EWMA-T-Squared Chart
#'
#' Estimates the in-control Average Run Length for an EWMA-smoothed
#' T-squared chart via Monte Carlo simulation.
#'
#' @param chart An object of class \code{spm.chart} from \code{\link{spm.phase1}}.
#' @param lambda EWMA smoothing parameter (default 0.2).
#' @param ucl Upper control limit for the EWMA chart. If \code{NULL}, uses
#'   \code{chart$t2.ucl}.
#' @param n.sim Number of simulations (default 10000).
#' @param max.rl Maximum run length per simulation (default 5000).
#' @param seed Random seed (default 42).
#'
#' @return An object of class \code{spm.arl} with the same structure as
#'   \code{\link{spm.arl}}.
#'
#' @seealso \code{\link{spm.arl}} for standard ARL, \code{\link{spm.ewma}}
#'   for EWMA monitoring
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 50; m <- 30
#' argvals <- seq(0, 1, length.out = m)
#' X <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(X, argvals = argvals)
#' chart <- spm.phase1(fd, ncomp = 3)
#'
#' arl_ewma <- spm.arl.ewma(chart, lambda = 0.2, n.sim = 1000)
#' arl_ewma
#' }
#'
#' @export
spm.arl.ewma <- function(chart, lambda = 0.2, ucl = NULL,
                           n.sim = 10000, max.rl = 5000, seed = 42) {
  if (!inherits(chart, "spm.chart")) {
    stop("chart must be of class 'spm.chart'")
  }

  if (is.null(ucl)) {
    ucl <- chart$t2.ucl
  }

  result <- .Call("wrap__spm_arl0_ewma_t2_rust",
                  as.numeric(chart$eigenvalues),
                  as.numeric(ucl), as.numeric(lambda),
                  as.integer(n.sim), as.integer(max.rl),
                  as.integer(seed))

  if (is.null(result)) {
    stop("spm.arl.ewma failed: check chart and parameters")
  }

  structure(
    list(
      arl = result$arl,
      std.dev = result$std_dev,
      median.rl = result$median_rl,
      run.lengths = result$run_lengths,
      type = "ewma.t2",
      in.control = TRUE,
      call = match.call()
    ),
    class = "spm.arl"
  )
}

#' CUSUM Monitoring for Functional Data
#'
#' Applies Cumulative Sum (CUSUM) monitoring to sequential functional data.
#' CUSUM accumulates evidence of a mean shift over time, providing faster
#' detection of small persistent shifts than Shewhart-type charts.
#'
#' @param chart An object of class \code{spm.chart} from \code{\link{spm.phase1}}.
#' @param newdata An object of class \code{fdata} with sequential observations.
#' @param k CUSUM allowance (slack) parameter (default 0.5). Controls the
#'   size of shift the CUSUM is tuned to detect.
#' @param h CUSUM decision interval (default 5.0). Larger values give fewer
#'   false alarms but slower detection.
#' @param ncomp Number of components for CUSUM (default: same as chart).
#' @param alpha Significance level for SPE limit (default 0.05).
#' @param restart Logical; if \code{TRUE}, the CUSUM resets to zero after each
#'   alarm (default \code{FALSE}).
#'
#' @return An object of class \code{spm.cusum} with components:
#' \describe{
#'   \item{cusum.statistic}{Numeric vector of CUSUM values}
#'   \item{ucl}{Decision interval (h)}
#'   \item{alarm}{Logical vector: TRUE where CUSUM exceeds h}
#'   \item{scores}{FPC score matrix}
#'   \item{spe}{SPE values}
#'   \item{spe.limit}{SPE control limit}
#'   \item{spe.alarm}{Logical: TRUE where SPE exceeds limit}
#' }
#'
#' @seealso \code{\link{spm.monitor}} for Shewhart-type monitoring,
#'   \code{\link{spm.mewma}} for MEWMA monitoring
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 50; m <- 30
#' argvals <- seq(0, 1, length.out = m)
#' X <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(X, argvals = argvals)
#' chart <- spm.phase1(fd, ncomp = 3)
#'
#' # Monitor with a small persistent shift
#' X_new <- matrix(rnorm(30 * m) + 0.5, 30, m)
#' fd_new <- fdata(X_new, argvals = argvals)
#' cusum <- spm.cusum(chart, fd_new, k = 0.5, h = 5.0)
#' cusum
#' }
#'
#' @export
spm.cusum <- function(chart, newdata, k = 0.5, h = 5.0,
                       ncomp = NULL, alpha = 0.05, restart = FALSE) {
  if (!inherits(chart, "spm.chart")) {
    stop("chart must be of class 'spm.chart'")
  }
  if (!inherits(newdata, "fdata")) {
    stop("newdata must be of class 'fdata'")
  }

  if (is.null(ncomp)) {
    ncomp <- chart$ncomp
  }

  r <- chart$.rust
  argvals <- as.numeric(chart$fdataobj$argvals)

  rust_fn <- if (restart) {
    "wrap__spm_cusum_monitor_restart_rust"
  } else {
    "wrap__spm_cusum_monitor_rust"
  }

  result <- .Call(rust_fn,
                  newdata$data, argvals,
                  r$rotation, r$mean, r$singular_values,
                  r$centered, r$eigenvalues,
                  r$t2_ucl, r$t2_alpha, r$t2_description,
                  r$spe_ucl, r$spe_alpha, r$spe_description,
                  as.integer(r$ncomp), r$config_alpha,
                  r$config_tuning_fraction, as.integer(r$config_seed),
                  as.numeric(k), as.numeric(h),
                  as.integer(ncomp), as.numeric(alpha),
                  FALSE)

  if (is.null(result)) {
    stop("spm.cusum failed: check data dimensions and parameters")
  }

  structure(
    list(
      cusum.statistic = result$cusum_statistic,
      ucl = result$ucl,
      alarm = as.logical(result$alarm),
      scores = result$scores,
      spe = result$spe,
      spe.limit = result$spe_limit,
      spe.alarm = as.logical(result$spe_alarm),
      k = k,
      h = h,
      restart = restart,
      call = match.call()
    ),
    class = "spm.cusum"
  )
}

#' @export
print.spm.cusum <- function(x, ...) {
  n <- length(x$cusum.statistic)
  n_alarm <- sum(x$alarm)
  n_spe <- sum(x$spe.alarm)

  cat("SPM CUSUM Monitoring\n")
  cat("  Observations:", n, "\n")
  cat("  k:", format(x$k, digits = 3), "\n")
  cat("  h (UCL):", format(x$ucl, digits = 4), "\n")
  cat("  Restart:", x$restart, "\n")
  cat("  CUSUM alarms:", n_alarm, "of", n,
      paste0("(", format(100 * n_alarm / n, digits = 3), "%)"), "\n")
  cat("  SPE alarms:", n_spe, "of", n,
      paste0("(", format(100 * n_spe / n, digits = 3), "%)"), "\n")
  invisible(x)
}

#' Plot CUSUM monitoring results
#'
#' Displays CUSUM statistic and SPE over time with control limits.
#'
#' @param x An object of class \code{spm.cusum}.
#' @param ... Additional arguments (ignored).
#' @export
plot.spm.cusum <- function(x, ...) {
  n <- length(x$cusum.statistic)
  df <- data.frame(
    obs = rep(seq_len(n), 2),
    value = c(x$cusum.statistic, x$spe),
    panel = factor(
      rep(c("CUSUM Statistic", "Squared Prediction Error"), each = n),
      levels = c("CUSUM Statistic", "Squared Prediction Error")
    ),
    ucl = rep(c(x$ucl, x$spe.limit), each = n),
    alarm = c(x$alarm, x$spe.alarm)
  )

  ucl_df <- data.frame(
    panel = factor(
      c("CUSUM Statistic", "Squared Prediction Error"),
      levels = c("CUSUM Statistic", "Squared Prediction Error")
    ),
    ucl = c(x$ucl, x$spe.limit),
    label = paste("UCL =", format(c(x$ucl, x$spe.limit), digits = 4))
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$obs, y = .data$value)) +
    ggplot2::geom_line(ggplot2::aes(color = .data$alarm), linewidth = 0.6) +
    ggplot2::geom_point(
      data = df[df$alarm, , drop = FALSE],
      ggplot2::aes(x = .data$obs, y = .data$value),
      color = "#D55E00", size = 1.5
    ) +
    ggplot2::geom_hline(
      data = ucl_df,
      ggplot2::aes(yintercept = .data$ucl),
      linetype = "dashed", color = "#D55E00", linewidth = 0.7
    ) +
    ggplot2::geom_text(
      data = ucl_df,
      ggplot2::aes(x = -Inf, y = .data$ucl, label = .data$label),
      hjust = -0.05, vjust = -0.5, size = 3, color = "#D55E00"
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "#4A90D9", "TRUE" = "#D55E00"),
      guide = "none"
    ) +
    ggplot2::facet_wrap(~ .data$panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = "CUSUM Monitoring",
      x = "Observation", y = "Statistic"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))

  p
}

# =============================================================================
# Phase 3: MEWMA, AMEWMA, Iterative Phase I
# =============================================================================

#' MEWMA Monitoring for Functional Data
#'
#' Applies Multivariate Exponentially Weighted Moving Average (MEWMA) monitoring
#' to sequential functional data. MEWMA smooths the FPC score vectors over time,
#' then computes a Hotelling-type statistic on the smoothed scores.
#'
#' @param chart An object of class \code{spm.chart} from \code{\link{spm.phase1}}.
#' @param newdata An object of class \code{fdata} with sequential observations.
#' @param lambda MEWMA smoothing parameter in (0, 1] (default 0.2).
#' @param ncomp Number of components (default: same as chart).
#' @param alpha Significance level (default 0.05).
#'
#' @return An object of class \code{spm.mewma} with components:
#' \describe{
#'   \item{smoothed.scores}{Matrix of MEWMA-smoothed score vectors}
#'   \item{mewma.statistic}{Numeric vector of MEWMA statistics}
#'   \item{ucl}{Control limit}
#'   \item{alarm}{Logical vector: TRUE where MEWMA exceeds UCL}
#'   \item{spe}{SPE values}
#'   \item{spe.limit}{SPE control limit}
#'   \item{spe.alarm}{Logical: TRUE where SPE exceeds limit}
#' }
#'
#' @seealso \code{\link{spm.amewma}} for adaptive MEWMA,
#'   \code{\link{spm.ewma}} for univariate EWMA
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 50; m <- 30
#' argvals <- seq(0, 1, length.out = m)
#' X <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(X, argvals = argvals)
#' chart <- spm.phase1(fd, ncomp = 3)
#'
#' X_new <- matrix(rnorm(20 * m) + 0.5, 20, m)
#' fd_new <- fdata(X_new, argvals = argvals)
#' mewma <- spm.mewma(chart, fd_new, lambda = 0.2)
#' mewma
#' }
#'
#' @export
spm.mewma <- function(chart, newdata, lambda = 0.2, ncomp = NULL,
                       alpha = 0.05) {
  if (!inherits(chart, "spm.chart")) {
    stop("chart must be of class 'spm.chart'")
  }
  if (!inherits(newdata, "fdata")) {
    stop("newdata must be of class 'fdata'")
  }

  if (is.null(ncomp)) {
    ncomp <- chart$ncomp
  }

  r <- chart$.rust
  argvals <- as.numeric(chart$fdataobj$argvals)

  result <- .Call("wrap__spm_mewma_monitor_rust",
                  newdata$data, argvals,
                  r$rotation, r$mean, r$singular_values,
                  r$centered, r$eigenvalues,
                  r$t2_ucl, r$t2_alpha, r$t2_description,
                  r$spe_ucl, r$spe_alpha, r$spe_description,
                  as.integer(r$ncomp), r$config_alpha,
                  r$config_tuning_fraction, as.integer(r$config_seed),
                  as.numeric(lambda), as.integer(ncomp),
                  as.numeric(alpha), TRUE)

  if (is.null(result)) {
    stop("spm.mewma failed: check data dimensions and parameters")
  }

  structure(
    list(
      smoothed.scores = result$smoothed_scores,
      mewma.statistic = result$mewma_statistic,
      ucl = result$ucl,
      alarm = as.logical(result$alarm),
      spe = result$spe,
      spe.limit = result$spe_limit,
      spe.alarm = as.logical(result$spe_alarm),
      lambda = lambda,
      call = match.call()
    ),
    class = "spm.mewma"
  )
}

#' @export
print.spm.mewma <- function(x, ...) {
  n <- length(x$mewma.statistic)
  n_alarm <- sum(x$alarm)
  n_spe <- sum(x$spe.alarm)

  cat("SPM MEWMA Monitoring\n")
  cat("  Observations:", n, "\n")
  cat("  Lambda:", format(x$lambda, digits = 3), "\n")
  cat("  UCL:", format(x$ucl, digits = 4), "\n")
  cat("  MEWMA alarms:", n_alarm, "of", n,
      paste0("(", format(100 * n_alarm / n, digits = 3), "%)"), "\n")
  cat("  SPE alarms:", n_spe, "of", n,
      paste0("(", format(100 * n_spe / n, digits = 3), "%)"), "\n")
  invisible(x)
}

#' Plot MEWMA monitoring results
#'
#' Displays MEWMA statistic and SPE over time with control limits.
#'
#' @param x An object of class \code{spm.mewma}.
#' @param ... Additional arguments (ignored).
#' @export
plot.spm.mewma <- function(x, ...) {
  n <- length(x$mewma.statistic)
  df <- data.frame(
    obs = rep(seq_len(n), 2),
    value = c(x$mewma.statistic, x$spe),
    panel = factor(
      rep(c("MEWMA Statistic", "Squared Prediction Error"), each = n),
      levels = c("MEWMA Statistic", "Squared Prediction Error")
    ),
    ucl = rep(c(x$ucl, x$spe.limit), each = n),
    alarm = c(x$alarm, x$spe.alarm)
  )

  ucl_df <- data.frame(
    panel = factor(
      c("MEWMA Statistic", "Squared Prediction Error"),
      levels = c("MEWMA Statistic", "Squared Prediction Error")
    ),
    ucl = c(x$ucl, x$spe.limit),
    label = paste("UCL =", format(c(x$ucl, x$spe.limit), digits = 4))
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$obs, y = .data$value)) +
    ggplot2::geom_line(ggplot2::aes(color = .data$alarm), linewidth = 0.6) +
    ggplot2::geom_point(
      data = df[df$alarm, , drop = FALSE],
      ggplot2::aes(x = .data$obs, y = .data$value),
      color = "#D55E00", size = 1.5
    ) +
    ggplot2::geom_hline(
      data = ucl_df,
      ggplot2::aes(yintercept = .data$ucl),
      linetype = "dashed", color = "#D55E00", linewidth = 0.7
    ) +
    ggplot2::geom_text(
      data = ucl_df,
      ggplot2::aes(x = -Inf, y = .data$ucl, label = .data$label),
      hjust = -0.05, vjust = -0.5, size = 3, color = "#D55E00"
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "#4A90D9", "TRUE" = "#D55E00"),
      guide = "none"
    ) +
    ggplot2::facet_wrap(~ .data$panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = "MEWMA Monitoring",
      x = "Observation", y = "Statistic"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))

  p
}

#' Adaptive MEWMA (AMEWMA) Monitoring for Functional Data
#'
#' Applies Adaptive Multivariate EWMA monitoring where the smoothing parameter
#' adapts over time based on the magnitude of recent deviations. This provides
#' good detection performance for both small and large shifts.
#'
#' @param chart An object of class \code{spm.chart} from \code{\link{spm.phase1}}.
#' @param newdata An object of class \code{fdata} with sequential observations.
#' @param lambda.min Minimum smoothing parameter (default 0.05).
#' @param lambda.max Maximum smoothing parameter (default 0.95).
#' @param lambda.init Initial smoothing parameter (default 0.2).
#' @param eta Learning rate for lambda adaptation (default 0.1).
#' @param ncomp Number of components (default: same as chart).
#' @param alpha Significance level (default 0.05).
#'
#' @return An object of class \code{spm.amewma} with components:
#' \describe{
#'   \item{smoothed.scores}{Matrix of adaptively smoothed score vectors}
#'   \item{t2.statistic}{Numeric vector of T-squared statistics on smoothed scores}
#'   \item{lambda.t}{Numeric vector of adaptive lambda values over time}
#'   \item{ucl}{Control limit}
#'   \item{alarm}{Logical vector: TRUE where statistic exceeds UCL}
#' }
#'
#' @seealso \code{\link{spm.mewma}} for fixed-lambda MEWMA,
#'   \code{\link{spm.ewma}} for univariate EWMA
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 50; m <- 30
#' argvals <- seq(0, 1, length.out = m)
#' X <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(X, argvals = argvals)
#' chart <- spm.phase1(fd, ncomp = 3)
#'
#' X_new <- matrix(rnorm(30 * m) + 0.5, 30, m)
#' fd_new <- fdata(X_new, argvals = argvals)
#' amewma <- spm.amewma(chart, fd_new)
#' amewma
#' }
#'
#' @export
spm.amewma <- function(chart, newdata, lambda.min = 0.05, lambda.max = 0.95,
                        lambda.init = 0.2, eta = 0.1,
                        ncomp = NULL, alpha = 0.05) {
  if (!inherits(chart, "spm.chart")) {
    stop("chart must be of class 'spm.chart'")
  }
  if (!inherits(newdata, "fdata")) {
    stop("newdata must be of class 'fdata'")
  }

  if (is.null(ncomp)) {
    ncomp <- chart$ncomp
  }

  r <- chart$.rust
  argvals <- as.numeric(chart$fdataobj$argvals)

  result <- .Call("wrap__spm_amewma_monitor_rust",
                  newdata$data, argvals,
                  r$rotation, r$mean, r$singular_values,
                  r$centered, r$eigenvalues,
                  r$t2_ucl, r$t2_alpha, r$t2_description,
                  r$spe_ucl, r$spe_alpha, r$spe_description,
                  as.integer(r$ncomp), r$config_alpha,
                  r$config_tuning_fraction, as.integer(r$config_seed),
                  as.numeric(lambda.min), as.numeric(lambda.max),
                  as.numeric(lambda.init), as.numeric(eta),
                  as.integer(ncomp), as.numeric(alpha))

  if (is.null(result)) {
    stop("spm.amewma failed: check data dimensions and parameters")
  }

  structure(
    list(
      smoothed.scores = result$smoothed_scores,
      t2.statistic = result$t2_statistic,
      lambda.t = result$lambda_t,
      ucl = result$ucl,
      alarm = as.logical(result$alarm),
      lambda.min = lambda.min,
      lambda.max = lambda.max,
      lambda.init = lambda.init,
      eta = eta,
      call = match.call()
    ),
    class = "spm.amewma"
  )
}

#' @export
print.spm.amewma <- function(x, ...) {
  n <- length(x$t2.statistic)
  n_alarm <- sum(x$alarm)

  cat("SPM Adaptive MEWMA (AMEWMA) Monitoring\n")
  cat("  Observations:", n, "\n")
  cat("  Lambda range: [", format(x$lambda.min, digits = 3), ", ",
      format(x$lambda.max, digits = 3), "]\n", sep = "")
  cat("  Lambda init:", format(x$lambda.init, digits = 3), "\n")
  cat("  Eta:", format(x$eta, digits = 3), "\n")
  cat("  UCL:", format(x$ucl, digits = 4), "\n")
  cat("  Alarms:", n_alarm, "of", n,
      paste0("(", format(100 * n_alarm / n, digits = 3), "%)"), "\n")
  invisible(x)
}

#' Plot AMEWMA monitoring results
#'
#' Displays AMEWMA T-squared statistic and adaptive lambda over time in
#' a two-panel plot.
#'
#' @param x An object of class \code{spm.amewma}.
#' @param ... Additional arguments (ignored).
#' @export
plot.spm.amewma <- function(x, ...) {
  n <- length(x$t2.statistic)
  obs <- seq_len(n)

  # Top panel: T2 statistic
  df_t2 <- data.frame(
    obs = obs,
    value = x$t2.statistic,
    alarm = x$alarm
  )

  p1 <- ggplot2::ggplot(df_t2, ggplot2::aes(x = .data$obs, y = .data$value)) +
    ggplot2::geom_line(ggplot2::aes(color = .data$alarm), linewidth = 0.6) +
    ggplot2::geom_point(
      data = df_t2[df_t2$alarm, , drop = FALSE],
      ggplot2::aes(x = .data$obs, y = .data$value),
      color = "#D55E00", size = 1.5
    ) +
    ggplot2::geom_hline(yintercept = x$ucl, linetype = "dashed",
                        color = "#D55E00", linewidth = 0.7) +
    ggplot2::annotate("text", x = -Inf, y = x$ucl,
                      label = paste("UCL =", format(x$ucl, digits = 4)),
                      hjust = -0.05, vjust = -0.5, size = 3,
                      color = "#D55E00") +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "#4A90D9", "TRUE" = "#D55E00"),
      guide = "none"
    ) +
    ggplot2::labs(title = "AMEWMA Monitoring",
                  x = NULL, y = "T\u00b2 Statistic") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())

  # Bottom panel: Adaptive lambda
  df_lambda <- data.frame(
    obs = obs,
    lambda = x$lambda.t
  )

  p2 <- ggplot2::ggplot(df_lambda, ggplot2::aes(x = .data$obs, y = .data$lambda)) +
    ggplot2::geom_line(color = "#4A90D9", linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = x$lambda.init, linetype = "dotted",
                        color = "grey50") +
    ggplot2::labs(x = "Observation", y = "Adaptive \u03bb") +
    ggplot2::ylim(x$lambda.min, x$lambda.max) +
    ggplot2::theme_minimal()

  p1 + p2 + patchwork::plot_layout(ncol = 1, heights = c(2, 1))
}

#' Iterative Phase I Chart Construction
#'
#' Builds a Phase I control chart iteratively by removing observations that
#' exceed control limits, then re-fitting. This produces cleaner reference
#' limits when the training data may contain outliers.
#'
#' @param fdataobj An object of class \code{fdata} with in-control data.
#' @param ncomp Number of principal components (default 5).
#' @param alpha Significance level (default 0.05).
#' @param tuning.fraction Fraction for FPCA tuning (default 0.5).
#' @param seed Random seed (default 42).
#' @param max.iterations Maximum number of iterative removal rounds (default 10).
#' @param max.removal.fraction Maximum fraction of data that may be removed
#'   (default 0.3).
#'
#' @return An object of class \code{spm.chart} (compatible with \code{spm.monitor},
#'   \code{spm.cusum}, etc.) with additional fields:
#'   \item{n.iterations}{Number of iterations performed}
#'   \item{removed.indices}{Integer vector of removed observation indices (1-indexed)}
#'   \item{n.remaining}{Number of observations remaining}
#'   \item{removal.history}{List of integer vectors; removed indices per iteration}
#'   \item{removal.rates}{Numeric vector of removal rates per iteration}
#'
#' @seealso \code{\link{spm.phase1}} for non-iterative Phase I
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 60; m <- 30
#' argvals <- seq(0, 1, length.out = m)
#' # Add some outliers
#' X <- matrix(rnorm(n * m), n, m)
#' X[1:3, ] <- X[1:3, ] + 5
#' fd <- fdata(X, argvals = argvals)
#'
#' chart <- spm.phase1.iterative(fd, ncomp = 3, max.iterations = 5)
#' chart
#' chart$removed.indices
#' }
#'
#' @export
spm.phase1.iterative <- function(fdataobj, ncomp = 5, alpha = 0.05,
                                  tuning.fraction = 0.5, seed = 42,
                                  max.iterations = 10,
                                  max.removal.fraction = 0.3) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  n <- nrow(fdataobj$data)
  if (n < 4) {
    stop("At least 4 observations are required")
  }

  argvals <- as.numeric(fdataobj$argvals)

  result <- .Call("wrap__spm_phase1_iterative_rust",
                  fdataobj$data, argvals,
                  as.integer(ncomp), as.numeric(alpha),
                  as.numeric(tuning.fraction), as.integer(seed),
                  as.integer(max.iterations),
                  as.numeric(max.removal.fraction))

  if (is.null(result)) {
    stop("spm.phase1.iterative failed: check data and parameters")
  }

  structure(
    list(
      eigenvalues = result$eigenvalues,
      t2.phase1 = result$t2_phase1,
      spe.phase1 = result$spe_phase1,
      t2.ucl = result$t2_ucl,
      spe.ucl = result$spe_ucl,
      ncomp = result$ncomp,
      alpha = alpha,
      t2.description = result$t2_description,
      spe.description = result$spe_description,
      fdataobj = fdataobj,
      n.iterations = result$n_iterations,
      removed.indices = result$removed_indices,
      n.remaining = result$n_remaining,
      removal.history = result$removal_history,
      removal.rates = result$removal_rates,
      call = match.call(),
      .rust = list(
        rotation = result$rotation,
        scores = result$scores,
        mean = result$mean,
        singular_values = result$singular_values,
        centered = result$centered,
        eigenvalues = result$eigenvalues,
        t2_ucl = result$t2_ucl,
        t2_alpha = result$t2_alpha,
        t2_description = result$t2_description,
        spe_ucl = result$spe_ucl,
        spe_alpha = result$spe_alpha,
        spe_description = result$spe_description,
        ncomp = result$ncomp,
        config_alpha = result$config_alpha,
        config_tuning_fraction = result$config_tuning_fraction,
        config_seed = result$config_seed
      )
    ),
    class = "spm.chart"
  )
}

# =============================================================================
# Phase 4: Profile monitoring and partial-domain monitoring
# =============================================================================

#' Profile Monitoring Phase I
#'
#' Builds a profile monitoring chart for functional data with scalar covariates.
#' Uses function-on-scalar regression (FOSR) to model the relationship, then
#' applies FPCA to the regression coefficients across sliding windows.
#'
#' @param fdataobj An object of class \code{fdata} (functional response).
#' @param predictors A matrix of scalar predictors (n x p).
#' @param ncomp Number of principal components for beta FPCA (default 3).
#' @param alpha Significance level (default 0.05).
#' @param fosr.lambda FOSR smoothing parameter (default 1e-4).
#' @param window.size Window size for sliding-window FOSR (default 20).
#' @param step.size Step size between windows (default 1).
#'
#' @return An object of class \code{spm.profile.chart} with components:
#' \describe{
#'   \item{eigenvalues}{Eigenvalues from beta FPCA}
#'   \item{t2.ucl}{T-squared control limit}
#'   \item{t2.description}{Description of the limit}
#'   \item{lag1.autocorrelation}{Lag-1 autocorrelation of T-squared values}
#'   \item{effective.n.windows}{Effective number of independent windows}
#'   \item{fdataobj}{Original fdata object}
#'   \item{predictors}{Original predictor matrix}
#'   \item{.rust}{Internal fields for Phase II monitoring}
#' }
#'
#' @seealso \code{\link{spm.profile.monitor}} for Phase II monitoring
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 80; m <- 20
#' argvals <- seq(0, 1, length.out = m)
#' X_pred <- cbind(rnorm(n), rnorm(n))
#' Y <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(Y, argvals = argvals)
#'
#' chart <- spm.profile.phase1(fd, X_pred, ncomp = 2, window.size = 15)
#' chart
#' }
#'
#' @export
spm.profile.phase1 <- function(fdataobj, predictors, ncomp = 3, alpha = 0.05,
                                 fosr.lambda = 1e-4, window.size = 20,
                                 step.size = 1) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  predictors <- as.matrix(predictors)
  n <- nrow(fdataobj$data)

  if (nrow(predictors) != n) {
    stop("Number of rows in predictors must equal number of curves")
  }

  argvals <- as.numeric(fdataobj$argvals)

  result <- .Call("wrap__spm_profile_phase1_rust",
                  fdataobj$data, predictors, argvals,
                  as.numeric(fosr.lambda), as.integer(ncomp),
                  as.numeric(alpha), as.integer(window.size),
                  as.integer(step.size))

  if (is.null(result)) {
    stop("spm.profile.phase1 failed: check data dimensions and parameters")
  }

  structure(
    list(
      eigenvalues = result$eigenvalues,
      t2.ucl = result$t2_ucl,
      t2.description = result$t2_description,
      lag1.autocorrelation = result$lag1_autocorrelation,
      effective.n.windows = result$effective_n_windows,
      ncomp = ncomp,
      alpha = alpha,
      fdataobj = fdataobj,
      predictors = predictors,
      call = match.call(),
      .rust = list(
        beta_rotation = result$beta_rotation,
        beta_scores = result$beta_scores,
        beta_mean = result$beta_mean,
        beta_singular_values = result$beta_singular_values,
        beta_centered = result$beta_centered,
        ref_fosr_intercept = result$ref_fosr_intercept,
        ref_fosr_beta = result$ref_fosr_beta,
        ref_fosr_lambda = result$ref_fosr_lambda,
        eigenvalues = result$eigenvalues,
        t2_ucl = result$t2_ucl,
        t2_alpha = result$t2_alpha,
        t2_description = result$t2_description,
        lag1_autocorrelation = result$lag1_autocorrelation,
        effective_n_windows = result$effective_n_windows,
        config_fosr_lambda = result$config_fosr_lambda,
        config_ncomp = result$config_ncomp,
        config_alpha = result$config_alpha,
        config_window_size = result$config_window_size,
        config_step_size = result$config_step_size
      )
    ),
    class = "spm.profile.chart"
  )
}

#' @export
print.spm.profile.chart <- function(x, ...) {
  cat("SPM Profile Monitoring Chart (Phase I)\n")
  cat("  Components:", x$ncomp, "\n")
  cat("  Alpha:", x$alpha, "\n")
  cat("  T2 UCL:", format(x$t2.ucl, digits = 4), "\n")
  cat("  Observations:", nrow(x$fdataobj$data), "\n")
  cat("  Predictors:", ncol(x$predictors), "\n")
  cat("  Lag-1 autocorrelation:", format(x$lag1.autocorrelation, digits = 4), "\n")
  cat("  Effective windows:", format(x$effective.n.windows, digits = 1), "\n")
  invisible(x)
}

#' Monitor New Data Against Profile Chart (Phase II)
#'
#' Projects new functional data and scalar predictors through a trained profile
#' monitoring chart. Computes sliding-window FOSR betas, projects onto the
#' reference beta FPCA, and evaluates T-squared statistics.
#'
#' @param chart An object of class \code{spm.profile.chart} from
#'   \code{\link{spm.profile.phase1}}.
#' @param newdata An object of class \code{fdata} with new functional response.
#' @param new.predictors A matrix of new scalar predictors.
#'
#' @return An object of class \code{spm.monitor} with components:
#' \describe{
#'   \item{betas}{Matrix of estimated beta coefficients per window}
#'   \item{t2}{T-squared values for each window}
#'   \item{t2.alarm}{Logical: TRUE where T-squared exceeds UCL}
#'   \item{beta.scores}{Beta FPC scores for each window}
#'   \item{t2.ucl}{T-squared control limit}
#' }
#'
#' @seealso \code{\link{spm.profile.phase1}} for building the chart
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 80; m <- 20
#' argvals <- seq(0, 1, length.out = m)
#' X_pred <- cbind(rnorm(n), rnorm(n))
#' Y <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(Y, argvals = argvals)
#' chart <- spm.profile.phase1(fd, X_pred, ncomp = 2, window.size = 15)
#'
#' # Monitor new data
#' X_new <- cbind(rnorm(30), rnorm(30))
#' Y_new <- matrix(rnorm(30 * m) + 1, 30, m)
#' fd_new <- fdata(Y_new, argvals = argvals)
#' mon <- spm.profile.monitor(chart, fd_new, X_new)
#' mon
#' }
#'
#' @export
spm.profile.monitor <- function(chart, newdata, new.predictors) {
  if (!inherits(chart, "spm.profile.chart")) {
    stop("chart must be of class 'spm.profile.chart'")
  }
  if (!inherits(newdata, "fdata")) {
    stop("newdata must be of class 'fdata'")
  }

  new.predictors <- as.matrix(new.predictors)
  r <- chart$.rust
  argvals <- as.numeric(chart$fdataobj$argvals)

  result <- .Call("wrap__spm_profile_monitor_rust",
                  newdata$data, new.predictors, argvals,
                  r$beta_rotation, r$beta_mean,
                  r$beta_singular_values, r$beta_centered,
                  r$ref_fosr_intercept, r$ref_fosr_beta,
                  r$ref_fosr_lambda,
                  r$eigenvalues,
                  r$t2_ucl, r$t2_alpha, r$t2_description,
                  r$lag1_autocorrelation, r$effective_n_windows,
                  as.numeric(r$config_fosr_lambda),
                  as.integer(r$config_ncomp), as.numeric(r$config_alpha),
                  as.integer(r$config_window_size),
                  as.integer(r$config_step_size))

  if (is.null(result)) {
    stop("spm.profile.monitor failed: check data dimensions")
  }

  structure(
    list(
      betas = result$betas,
      t2 = result$t2,
      t2.alarm = as.logical(result$t2_alarm),
      beta.scores = result$beta_scores,
      t2.ucl = chart$t2.ucl,
      call = match.call()
    ),
    class = "spm.monitor"
  )
}

#' Monitor a Partially-Observed Curve
#'
#' Monitors a single curve that is only partially observed on the domain.
#' The unobserved portion is completed using one of several methods before
#' projecting onto the Phase I FPCA basis.
#'
#' @param chart An object of class \code{spm.chart} from \code{\link{spm.phase1}}.
#' @param partial.data Numeric vector of observed values for one curve.
#' @param n.observed Integer; how many of the grid points have been observed
#'   (from the start of the domain).
#' @param completion Completion method: \code{"conditional"} (conditional
#'   expectation / BLUP), \code{"projection"} (partial projection), or
#'   \code{"zero_pad"}.
#' @param ncomp Number of components (default: same as chart).
#' @param alpha Significance level (default 0.05).
#'
#' @return A list with components:
#' \describe{
#'   \item{scores}{Estimated FPC scores}
#'   \item{t2}{T-squared statistic}
#'   \item{t2.alarm}{Logical: TRUE if T-squared exceeds UCL}
#'   \item{domain.fraction}{Fraction of domain observed}
#'   \item{completed.curve}{The completed curve (full domain), or NULL}
#' }
#'
#' @seealso \code{\link{spm.monitor}} for fully-observed curves,
#'   \code{\link{spm.monitor.partial.batch}} for multiple partial curves
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 50; m <- 30
#' argvals <- seq(0, 1, length.out = m)
#' X <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(X, argvals = argvals)
#' chart <- spm.phase1(fd, ncomp = 3)
#'
#' # Monitor a partially observed curve (first 20 of 30 grid points)
#' partial <- rnorm(30)
#' result <- spm.monitor.partial(chart, partial, n.observed = 20)
#' result$t2
#' result$domain.fraction
#' }
#'
#' @export
spm.monitor.partial <- function(chart, partial.data, n.observed,
                                 completion = c("conditional", "projection",
                                                "zero_pad"),
                                 ncomp = NULL, alpha = 0.05) {
  if (!inherits(chart, "spm.chart")) {
    stop("chart must be of class 'spm.chart'")
  }

  completion <- match.arg(completion)
  partial.data <- as.numeric(partial.data)
  n.observed <- as.integer(n.observed)

  if (is.null(ncomp)) {
    ncomp <- chart$ncomp
  }

  # Map R completion names to Rust
  completion_rust <- switch(completion,
    conditional = "conditional_expectation",
    projection = "partial_projection",
    zero_pad = "zero_pad"
  )

  r <- chart$.rust
  argvals <- as.numeric(chart$fdataobj$argvals)

  result <- .Call("wrap__spm_monitor_partial_rust",
                  partial.data, argvals, n.observed,
                  r$rotation, r$mean, r$singular_values,
                  r$centered, r$eigenvalues,
                  r$t2_ucl, r$t2_alpha, r$t2_description,
                  r$spe_ucl, r$spe_alpha, r$spe_description,
                  as.integer(r$ncomp), r$config_alpha,
                  r$config_tuning_fraction, as.integer(r$config_seed),
                  as.integer(ncomp), as.numeric(alpha),
                  completion_rust)

  if (is.null(result)) {
    stop("spm.monitor.partial failed: check input dimensions")
  }

  list(
    scores = result$scores,
    t2 = result$t2,
    t2.alarm = as.logical(result$t2_alarm),
    domain.fraction = result$domain_fraction,
    completed.curve = result$completed_curve
  )
}

#' Monitor a Batch of Partially-Observed Curves
#'
#' Monitors multiple partially-observed curves in a single call. Each curve
#' may have a different number of observed grid points.
#'
#' @param chart An object of class \code{spm.chart} from \code{\link{spm.phase1}}.
#' @param partial.data.list A list of numeric vectors, each containing the
#'   observed values for one curve.
#' @param n.observed.list An integer vector giving the number of observed
#'   grid points for each curve.
#' @param completion Completion method: \code{"conditional"}, \code{"projection"},
#'   or \code{"zero_pad"}.
#' @param ncomp Number of components (default: same as chart).
#' @param alpha Significance level (default 0.05).
#'
#' @return A list of results, one per curve. Each element is a list with:
#'   \item{scores}{Estimated FPC scores}
#'   \item{t2}{T-squared statistic}
#'   \item{t2.alarm}{Logical}
#'   \item{domain.fraction}{Fraction of domain observed}
#'   \item{completed.curve}{The completed curve, or NULL}
#'
#' @seealso \code{\link{spm.monitor.partial}} for a single curve
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 50; m <- 30
#' argvals <- seq(0, 1, length.out = m)
#' X <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(X, argvals = argvals)
#' chart <- spm.phase1(fd, ncomp = 3)
#'
#' # Three curves with different observed lengths
#' curves <- list(rnorm(30), rnorm(30), rnorm(30))
#' n_obs <- c(10L, 20L, 25L)
#' results <- spm.monitor.partial.batch(chart, curves, n_obs)
#' sapply(results, function(r) r$t2)
#' }
#'
#' @export
spm.monitor.partial.batch <- function(chart, partial.data.list, n.observed.list,
                                       completion = c("conditional", "projection",
                                                      "zero_pad"),
                                       ncomp = NULL, alpha = 0.05) {
  if (!inherits(chart, "spm.chart")) {
    stop("chart must be of class 'spm.chart'")
  }
  if (!is.list(partial.data.list)) {
    stop("partial.data.list must be a list of numeric vectors")
  }

  completion <- match.arg(completion)
  n.observed.list <- as.integer(n.observed.list)

  if (length(partial.data.list) != length(n.observed.list)) {
    stop("partial.data.list and n.observed.list must have the same length")
  }

  if (is.null(ncomp)) {
    ncomp <- chart$ncomp
  }

  completion_rust <- switch(completion,
    conditional = "conditional_expectation",
    projection = "partial_projection",
    zero_pad = "zero_pad"
  )

  # Ensure all elements are numeric
  partial.data.list <- lapply(partial.data.list, as.numeric)

  r <- chart$.rust
  argvals <- as.numeric(chart$fdataobj$argvals)

  result <- .Call("wrap__spm_monitor_partial_batch_rust",
                  partial.data.list, n.observed.list, argvals,
                  r$rotation, r$mean, r$singular_values,
                  r$centered, r$eigenvalues,
                  r$t2_ucl, r$t2_alpha, r$t2_description,
                  r$spe_ucl, r$spe_alpha, r$spe_description,
                  as.integer(r$ncomp), r$config_alpha,
                  r$config_tuning_fraction, as.integer(r$config_seed),
                  as.integer(ncomp), as.numeric(alpha),
                  completion_rust)

  if (is.null(result)) {
    stop("spm.monitor.partial.batch failed: check input dimensions")
  }

  # Convert each element's t2_alarm to logical
  lapply(result, function(res) {
    list(
      scores = res$scores,
      t2 = res$t2,
      t2.alarm = as.logical(res$t2_alarm),
      domain.fraction = res$domain_fraction,
      completed.curve = res$completed_curve
    )
  })
}

# =============================================================================
# Phase 5: Elastic SPM
# =============================================================================

#' Elastic SPM Phase I Chart
#'
#' Builds an SPM control chart that is invariant to phase (time-warping)
#' variation. First aligns the functional data to a Karcher mean using
#' elastic registration, then builds separate amplitude and (optionally)
#' phase charts.
#'
#' @param fdataobj An object of class \code{fdata} with in-control data.
#' @param ncomp Number of principal components for amplitude chart (default 5).
#' @param alpha Significance level (default 0.05).
#' @param tuning.fraction Fraction for FPCA tuning (default 0.5).
#' @param seed Random seed (default 42).
#' @param align.lambda Elastic alignment regularization (default 0.0).
#' @param monitor.phase Logical; if \code{TRUE}, also builds a phase chart
#'   for warping function monitoring (default \code{TRUE}).
#' @param warp.ncomp Number of components for the phase (warping) chart
#'   (default 3).
#'
#' @return An object of class \code{spm.elastic.chart} with components:
#' \describe{
#'   \item{karcher.mean}{The Karcher mean function}
#'   \item{mean.alignment.residual}{Residual from Karcher mean computation}
#'   \item{amplitude}{List with amplitude chart fields: eigenvalues, t2.phase1,
#'     spe.phase1, t2.ucl, spe.ucl, ncomp}
#'   \item{phase}{List with phase chart fields (if \code{monitor.phase = TRUE}),
#'     or NULL}
#'   \item{fdataobj}{Original fdata object}
#'   \item{.rust}{Internal fields for Phase II monitoring}
#' }
#'
#' @seealso \code{\link{spm.elastic.monitor}} for Phase II elastic monitoring,
#'   \code{\link{spm.phase1}} for non-elastic Phase I
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 40; m <- 30
#' argvals <- seq(0, 1, length.out = m)
#' X <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(X, argvals = argvals)
#'
#' chart <- spm.elastic.phase1(fd, ncomp = 3, warp.ncomp = 2)
#' chart
#' }
#'
#' @export
spm.elastic.phase1 <- function(fdataobj, ncomp = 5, alpha = 0.05,
                                 tuning.fraction = 0.5, seed = 42,
                                 align.lambda = 0.0, monitor.phase = TRUE,
                                 warp.ncomp = 3) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  n <- nrow(fdataobj$data)
  if (n < 4) {
    stop("At least 4 observations are required")
  }

  argvals <- as.numeric(fdataobj$argvals)

  result <- .Call("wrap__elastic_spm_phase1_rust",
                  fdataobj$data, argvals,
                  as.integer(ncomp), as.numeric(alpha),
                  as.numeric(tuning.fraction), as.integer(seed),
                  as.numeric(align.lambda), as.logical(monitor.phase),
                  as.integer(warp.ncomp))

  if (is.null(result)) {
    stop("spm.elastic.phase1 failed: check data and parameters")
  }

  # Extract phase chart from attribute if present
  phase_chart_raw <- attr(result, "phase_chart")
  phase_info <- NULL
  if (!is.null(phase_chart_raw)) {
    phase_info <- list(
      eigenvalues = phase_chart_raw$eigenvalues,
      t2.phase1 = phase_chart_raw$t2_phase1,
      spe.phase1 = phase_chart_raw$spe_phase1,
      t2.ucl = phase_chart_raw$t2_ucl,
      spe.ucl = phase_chart_raw$spe_ucl,
      ncomp = phase_chart_raw$ncomp
    )
  }

  structure(
    list(
      karcher.mean = result$karcher_mean,
      mean.alignment.residual = result$mean_alignment_residual,
      amplitude = list(
        eigenvalues = result$amp_eigenvalues,
        t2.phase1 = result$amp_t2_phase1,
        spe.phase1 = result$amp_spe_phase1,
        t2.ucl = result$amp_t2_ucl,
        spe.ucl = result$amp_spe_ucl,
        ncomp = result$amp_ncomp
      ),
      phase = phase_info,
      ncomp = ncomp,
      alpha = alpha,
      monitor.phase = monitor.phase,
      warp.ncomp = warp.ncomp,
      fdataobj = fdataobj,
      call = match.call(),
      .rust = list(
        karcher_mean = result$karcher_mean,
        align_lambda = result$config_align_lambda,
        monitor_phase = result$config_monitor_phase,
        warp_ncomp = result$config_warp_ncomp,
        # Amplitude chart reconstruction
        amp_rotation = result$amp_rotation,
        amp_scores = result$amp_scores,
        amp_mean = result$amp_mean,
        amp_singular_values = result$amp_singular_values,
        amp_centered = result$amp_centered,
        amp_eigenvalues = result$amp_eigenvalues,
        amp_t2_ucl = result$amp_t2_ucl,
        amp_t2_alpha = result$amp_t2_alpha,
        amp_t2_description = result$amp_t2_description,
        amp_spe_ucl = result$amp_spe_ucl,
        amp_spe_alpha = result$amp_spe_alpha,
        amp_spe_description = result$amp_spe_description,
        amp_ncomp = result$amp_ncomp,
        config_ncomp = result$config_ncomp,
        config_alpha = result$config_alpha,
        config_tuning_fraction = result$config_tuning_fraction,
        config_seed = result$config_seed,
        # Phase chart reconstruction (may be NULL)
        phase_chart = phase_chart_raw
      )
    ),
    class = "spm.elastic.chart"
  )
}

#' @export
print.spm.elastic.chart <- function(x, ...) {
  cat("Elastic SPM Control Chart (Phase I)\n")
  cat("  Components (amplitude):", x$amplitude$ncomp, "\n")
  cat("  Alpha:", x$alpha, "\n")
  cat("  Amplitude T2 UCL:", format(x$amplitude$t2.ucl, digits = 4), "\n")
  cat("  Amplitude SPE UCL:", format(x$amplitude$spe.ucl, digits = 4), "\n")
  if (!is.null(x$phase)) {
    cat("  Phase monitoring: enabled\n")
    cat("  Phase components:", x$phase$ncomp, "\n")
    cat("  Phase T2 UCL:", format(x$phase$t2.ucl, digits = 4), "\n")
    cat("  Phase SPE UCL:", format(x$phase$spe.ucl, digits = 4), "\n")
  } else {
    cat("  Phase monitoring: disabled\n")
  }
  cat("  Observations:", nrow(x$fdataobj$data), "\n")
  cat("  Grid points:", ncol(x$fdataobj$data), "\n")
  invisible(x)
}

#' Plot Elastic SPM Phase I Chart
#'
#' Displays amplitude and (optionally) phase control charts for elastic
#' SPM Phase I.
#'
#' @param x An object of class \code{spm.elastic.chart}.
#' @param ... Additional arguments (ignored).
#' @export
plot.spm.elastic.chart <- function(x, ...) {
  # Amplitude chart
  amp <- x$amplitude
  n_amp <- length(amp$t2.phase1)
  df_amp <- data.frame(
    obs = rep(seq_len(n_amp), 2),
    value = c(amp$t2.phase1, amp$spe.phase1),
    panel = factor(
      rep(c("Amplitude T\u00b2", "Amplitude SPE"), each = n_amp),
      levels = c("Amplitude T\u00b2", "Amplitude SPE")
    ),
    ucl = rep(c(amp$t2.ucl, amp$spe.ucl), each = n_amp)
  )
  df_amp$alarm <- df_amp$value > df_amp$ucl

  ucl_amp <- data.frame(
    panel = factor(
      c("Amplitude T\u00b2", "Amplitude SPE"),
      levels = c("Amplitude T\u00b2", "Amplitude SPE")
    ),
    ucl = c(amp$t2.ucl, amp$spe.ucl),
    label = paste("UCL =", format(c(amp$t2.ucl, amp$spe.ucl), digits = 4))
  )

  p_amp <- ggplot2::ggplot(df_amp, ggplot2::aes(x = .data$obs, y = .data$value)) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = .data$obs, y = 0, yend = .data$value,
                   color = .data$alarm),
      linewidth = 0.8
    ) +
    ggplot2::geom_hline(
      data = ucl_amp,
      ggplot2::aes(yintercept = .data$ucl),
      linetype = "dashed", color = "#D55E00", linewidth = 0.7
    ) +
    ggplot2::geom_text(
      data = ucl_amp,
      ggplot2::aes(x = -Inf, y = .data$ucl, label = .data$label),
      hjust = -0.05, vjust = -0.5, size = 3, color = "#D55E00"
    ) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "#4A90D9", "TRUE" = "#D55E00"),
      guide = "none"
    ) +
    ggplot2::facet_wrap(~ .data$panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(
      title = "Elastic SPM Phase I: Amplitude",
      x = "Observation", y = "Statistic"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))

  if (!is.null(x$phase)) {
    ph <- x$phase
    n_ph <- length(ph$t2.phase1)
    df_ph <- data.frame(
      obs = rep(seq_len(n_ph), 2),
      value = c(ph$t2.phase1, ph$spe.phase1),
      panel = factor(
        rep(c("Phase T\u00b2", "Phase SPE"), each = n_ph),
        levels = c("Phase T\u00b2", "Phase SPE")
      ),
      ucl = rep(c(ph$t2.ucl, ph$spe.ucl), each = n_ph)
    )
    df_ph$alarm <- df_ph$value > df_ph$ucl

    ucl_ph <- data.frame(
      panel = factor(
        c("Phase T\u00b2", "Phase SPE"),
        levels = c("Phase T\u00b2", "Phase SPE")
      ),
      ucl = c(ph$t2.ucl, ph$spe.ucl),
      label = paste("UCL =", format(c(ph$t2.ucl, ph$spe.ucl), digits = 4))
    )

    p_ph <- ggplot2::ggplot(df_ph, ggplot2::aes(x = .data$obs, y = .data$value)) +
      ggplot2::geom_segment(
        ggplot2::aes(xend = .data$obs, y = 0, yend = .data$value,
                     color = .data$alarm),
        linewidth = 0.8
      ) +
      ggplot2::geom_hline(
        data = ucl_ph,
        ggplot2::aes(yintercept = .data$ucl),
        linetype = "dashed", color = "#D55E00", linewidth = 0.7
      ) +
      ggplot2::geom_text(
        data = ucl_ph,
        ggplot2::aes(x = -Inf, y = .data$ucl, label = .data$label),
        hjust = -0.05, vjust = -0.5, size = 3, color = "#D55E00"
      ) +
      ggplot2::scale_color_manual(
        values = c("FALSE" = "#4A90D9", "TRUE" = "#D55E00"),
        guide = "none"
      ) +
      ggplot2::facet_wrap(~ .data$panel, ncol = 1, scales = "free_y") +
      ggplot2::labs(
        title = "Elastic SPM Phase I: Phase",
        x = "Observation", y = "Statistic"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))

    return(p_amp + p_ph + patchwork::plot_layout(ncol = 1))
  }

  p_amp
}

#' Elastic SPM Phase II Monitoring
#'
#' Monitors new functional data against an elastic SPM chart. New curves are
#' aligned to the Karcher mean, then amplitude (and optionally phase) are
#' monitored separately.
#'
#' @param chart An object of class \code{spm.elastic.chart} from
#'   \code{\link{spm.elastic.phase1}}.
#' @param newdata An object of class \code{fdata} with new observations.
#'
#' @return A list with components:
#' \describe{
#'   \item{amplitude}{List with t2, spe, t2.alarm, spe.alarm, scores}
#'   \item{phase}{List with t2, spe, t2.alarm, spe.alarm, scores (if
#'     \code{monitor.phase = TRUE}), or NULL}
#'   \item{aligned.data}{Matrix of aligned curves}
#'   \item{warping.functions}{Matrix of estimated warping functions}
#' }
#'
#' @seealso \code{\link{spm.elastic.phase1}} for building the chart
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 40; m <- 30
#' argvals <- seq(0, 1, length.out = m)
#' X <- matrix(rnorm(n * m), n, m)
#' fd <- fdata(X, argvals = argvals)
#' chart <- spm.elastic.phase1(fd, ncomp = 3, warp.ncomp = 2)
#'
#' X_new <- matrix(rnorm(10 * m) + 1, 10, m)
#' fd_new <- fdata(X_new, argvals = argvals)
#' mon <- spm.elastic.monitor(chart, fd_new)
#' mon$amplitude$t2.alarm
#' }
#'
#' @export
spm.elastic.monitor <- function(chart, newdata) {
  if (!inherits(chart, "spm.elastic.chart")) {
    stop("chart must be of class 'spm.elastic.chart'")
  }
  if (!inherits(newdata, "fdata")) {
    stop("newdata must be of class 'fdata'")
  }

  r <- chart$.rust
  argvals <- as.numeric(chart$fdataobj$argvals)

  monitor_phase <- as.logical(r$monitor_phase)

  # Phase chart fields (use actual values or placeholders)
  pc <- r$phase_chart
  if (!is.null(pc)) {
    phase_rotation <- pc$rotation
    phase_mean <- pc$mean
    phase_singular_values <- pc$singular_values
    phase_centered <- pc$centered
    phase_eigenvalues <- pc$eigenvalues
    phase_t2_ucl <- pc$t2_ucl
    phase_t2_alpha <- pc$t2_alpha
    phase_t2_description <- pc$t2_description
    phase_spe_ucl <- pc$spe_ucl
    phase_spe_alpha <- pc$spe_alpha
    phase_spe_description <- pc$spe_description
    phase_ncomp <- as.integer(pc$ncomp)
    phase_config_alpha <- r$config_alpha
    phase_config_tuning_fraction <- r$config_tuning_fraction
    phase_config_seed <- as.integer(r$config_seed)
  } else {
    # Placeholders for when phase monitoring is disabled
    phase_rotation <- NULL
    phase_mean <- NULL
    phase_singular_values <- NULL
    phase_centered <- NULL
    phase_eigenvalues <- NULL
    phase_t2_ucl <- 0.0
    phase_t2_alpha <- 0.05
    phase_t2_description <- ""
    phase_spe_ucl <- 0.0
    phase_spe_alpha <- 0.05
    phase_spe_description <- ""
    phase_ncomp <- 1L
    phase_config_alpha <- r$config_alpha
    phase_config_tuning_fraction <- r$config_tuning_fraction
    phase_config_seed <- as.integer(r$config_seed)
  }

  result <- .Call("wrap__elastic_spm_monitor_rust",
                  newdata$data, argvals,
                  as.numeric(r$karcher_mean),
                  as.numeric(r$align_lambda),
                  monitor_phase,
                  as.integer(r$warp_ncomp),
                  # Amplitude chart
                  r$amp_rotation, r$amp_mean,
                  r$amp_singular_values, r$amp_centered,
                  r$amp_eigenvalues,
                  r$amp_t2_ucl, r$amp_t2_alpha,
                  r$amp_t2_description,
                  r$amp_spe_ucl, r$amp_spe_alpha,
                  r$amp_spe_description,
                  as.integer(r$amp_ncomp),
                  r$config_alpha, r$config_tuning_fraction,
                  as.integer(r$config_seed),
                  # Phase chart
                  phase_rotation, phase_mean,
                  phase_singular_values, phase_centered,
                  phase_eigenvalues,
                  phase_t2_ucl, phase_t2_alpha,
                  phase_t2_description,
                  phase_spe_ucl, phase_spe_alpha,
                  phase_spe_description,
                  phase_ncomp, phase_config_alpha,
                  phase_config_tuning_fraction,
                  phase_config_seed)

  if (is.null(result)) {
    stop("spm.elastic.monitor failed: check data dimensions")
  }

  # Extract phase results from attribute if present
  phase_raw <- attr(result, "phase")
  phase_result <- NULL
  if (!is.null(phase_raw)) {
    phase_result <- list(
      t2 = phase_raw$t2,
      spe = phase_raw$spe,
      t2.alarm = as.logical(phase_raw$t2_alarm),
      spe.alarm = as.logical(phase_raw$spe_alarm),
      scores = phase_raw$scores
    )
  }

  list(
    amplitude = list(
      t2 = result$amp_t2,
      spe = result$amp_spe,
      t2.alarm = as.logical(result$amp_t2_alarm),
      spe.alarm = as.logical(result$amp_spe_alarm),
      scores = result$amp_scores
    ),
    phase = phase_result,
    aligned.data = result$aligned_data,
    warping.functions = result$warping_functions
  )
}
