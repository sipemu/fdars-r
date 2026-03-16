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
#'   \item{t2}{T-squared values for new observations}
#'   \item{spe}{SPE values for new observations}
#'   \item{t2.alarm}{Logical vector: TRUE where T-squared exceeds UCL}
#'   \item{spe.alarm}{Logical vector: TRUE where SPE exceeds UCL}
#'   \item{scores}{FPC score matrix for new observations}
#'   \item{t2.ucl}{T-squared control limit}
#'   \item{spe.ucl}{SPE control limit}
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
#'   \item{scores}{Score matrix (n x ncomp)}
#'   \item{eigenfunctions}{List of eigenfunction matrices, one per variable}
#'   \item{eigenvalues}{Eigenvalues (length ncomp)}
#'   \item{means}{List of mean functions, one per variable}
#'   \item{scales}{Per-variable standard deviations}
#'   \item{grid.sizes}{Grid sizes per variable}
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
      .rust = list(combined_rotation = result$combined_rotation)
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
#'   \item{smoothed.scores}{EWMA-smoothed score matrix}
#'   \item{t2}{T-squared values on smoothed scores}
#'   \item{spe}{SPE values (unsmoothed)}
#'   \item{t2.limit}{T-squared control limit for EWMA}
#'   \item{spe.limit}{SPE control limit}
#'   \item{t2.alarm}{Logical: TRUE where T-squared exceeds limit}
#'   \item{spe.alarm}{Logical: TRUE where SPE exceeds limit}
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
#'   \item{eigenvalues}{Eigenvalues from residual FPCA}
#'   \item{t2.ucl}{T-squared control limit}
#'   \item{spe.ucl}{SPE control limit}
#'   \item{ncomp}{Number of components used}
#'   \item{fdataobj}{Original fdata object}
#'   \item{predictors}{Original predictor matrix}
#'   \item{.rust}{Internal fields for Phase II monitoring}
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
#'   \item{t2}{T-squared values}
#'   \item{spe}{SPE values}
#'   \item{t2.alarm}{Logical: TRUE where T-squared exceeds UCL}
#'   \item{spe.alarm}{Logical: TRUE where SPE exceeds UCL}
#'   \item{residual.scores}{Residual FPC scores}
#'   \item{t2.ucl}{T-squared control limit}
#'   \item{spe.ucl}{SPE control limit}
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
