#' Conformal Prediction for Functional Linear Model
#'
#' Constructs prediction intervals using split conformal inference for
#' the FPC-based functional linear model.
#'
#' @param fdataobj An object of class 'fdata' (training data).
#' @param y Response vector (training).
#' @param newdata An object of class 'fdata' (test data).
#' @param scalar.train Optional scalar covariates for training.
#' @param scalar.test Optional scalar covariates for test.
#' @param ncomp Number of FPC components (default 3).
#' @param cal.fraction Fraction of data for calibration (default 0.25).
#' @param alpha Miscoverage level (default 0.1 for 90 percent intervals).
#' @param seed Random seed.
#'
#' @return A list with components:
#' \describe{
#'   \item{predictions}{Point predictions for test data}
#'   \item{lower}{Lower bounds of prediction intervals}
#'   \item{upper}{Upper bounds of prediction intervals}
#'   \item{residual.quantile}{Calibration residual quantile}
#'   \item{coverage}{Empirical coverage on calibration set}
#' }
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' cp <- conformal.fregre.lm(fd[1:40, ], y[1:40], fd[41:50, ], ncomp = 3)
#' }
#'
#' @export
conformal.fregre.lm <- function(fdataobj, y, newdata,
                                 scalar.train = NULL, scalar.test = NULL,
                                 ncomp = 3, cal.fraction = 0.25,
                                 alpha = 0.1, seed = NULL) {
  if (!inherits(fdataobj, "fdata")) stop("fdataobj must be of class 'fdata'")
  if (!inherits(newdata, "fdata")) stop("newdata must be of class 'fdata'")
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  st <- if (!is.null(scalar.train)) as.matrix(scalar.train) else NULL
  ste <- if (!is.null(scalar.test)) as.matrix(scalar.test) else NULL

  result <- .Call("wrap__conformal_fregre_lm_rust",
                  fdataobj$data, as.numeric(y), newdata$data,
                  st, ste, as.integer(ncomp),
                  as.numeric(cal.fraction), as.numeric(alpha),
                  as.integer(seed))

  if (is.null(result)) stop("conformal.fregre.lm failed")
  result
}

#' Conformal Prediction for Nonparametric Functional Regression
#'
#' @param fdataobj An object of class 'fdata' (training data).
#' @param y Response vector (training).
#' @param newdata An object of class 'fdata' (test data).
#' @param scalar.train Optional scalar covariates for training.
#' @param scalar.test Optional scalar covariates for test.
#' @param h.func Functional bandwidth (0 for auto).
#' @param h.scalar Scalar bandwidth (0 for auto).
#' @param cal.fraction Fraction of data for calibration (default 0.25).
#' @param alpha Miscoverage level (default 0.1).
#' @param seed Random seed.
#'
#' @return Same as \code{\link{conformal.fregre.lm}}.
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' cp <- conformal.fregre.np(fd[1:40, ], y[1:40], fd[41:50, ])
#' }
#'
#' @export
conformal.fregre.np <- function(fdataobj, y, newdata,
                                 scalar.train = NULL, scalar.test = NULL,
                                 h.func = 0, h.scalar = 0,
                                 cal.fraction = 0.25, alpha = 0.1,
                                 seed = NULL) {
  if (!inherits(fdataobj, "fdata")) stop("fdataobj must be of class 'fdata'")
  if (!inherits(newdata, "fdata")) stop("newdata must be of class 'fdata'")
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  st <- if (!is.null(scalar.train)) as.matrix(scalar.train) else NULL
  ste <- if (!is.null(scalar.test)) as.matrix(scalar.test) else NULL

  result <- .Call("wrap__conformal_fregre_np_rust",
                  fdataobj$data, as.numeric(y), newdata$data,
                  as.numeric(fdataobj$argvals), st, ste,
                  as.numeric(h.func), as.numeric(h.scalar),
                  as.numeric(cal.fraction), as.numeric(alpha),
                  as.integer(seed))

  if (is.null(result)) stop("conformal.fregre.np failed")
  result
}

#' Conformal Prediction for Elastic Regression
#'
#' @param fdataobj An object of class 'fdata' (training data).
#' @param y Response vector (training).
#' @param newdata An object of class 'fdata' (test data).
#' @param ncomp.beta Number of basis functions for beta (default 10).
#' @param lambda Regularization (default 0).
#' @param cal.fraction Calibration fraction (default 0.25).
#' @param alpha Miscoverage level (default 0.1).
#' @param seed Random seed.
#'
#' @return Same as \code{\link{conformal.fregre.lm}}.
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' cp <- conformal.elastic.regression(fd[1:40, ], y[1:40], fd[41:50, ])
#' }
#'
#' @export
conformal.elastic.regression <- function(fdataobj, y, newdata,
                                          ncomp.beta = 10, lambda = 0,
                                          cal.fraction = 0.25, alpha = 0.1,
                                          seed = NULL) {
  if (!inherits(fdataobj, "fdata")) stop("fdataobj must be of class 'fdata'")
  if (!inherits(newdata, "fdata")) stop("newdata must be of class 'fdata'")
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  result <- .Call("wrap__conformal_elastic_regression_rust",
                  fdataobj$data, as.numeric(y), newdata$data,
                  as.numeric(fdataobj$argvals),
                  as.integer(ncomp.beta), as.numeric(lambda),
                  as.numeric(cal.fraction), as.numeric(alpha),
                  as.integer(seed))

  if (is.null(result)) stop("conformal.elastic.regression failed")
  result
}

#' Conformal Prediction for Elastic PCR
#'
#' @param fdataobj An object of class 'fdata' (training data).
#' @param y Response vector (training).
#' @param newdata An object of class 'fdata' (test data).
#' @param ncomp Number of FPC components (default 3).
#' @param pca.method PCA method: "vertical", "horizontal", or "joint".
#' @param lambda Regularization (default 0).
#' @param cal.fraction Calibration fraction (default 0.25).
#' @param alpha Miscoverage level (default 0.1).
#' @param seed Random seed.
#'
#' @return Same as \code{\link{conformal.fregre.lm}}.
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' cp <- conformal.elastic.pcr(fd[1:40, ], y[1:40], fd[41:50, ])
#' }
#'
#' @export
conformal.elastic.pcr <- function(fdataobj, y, newdata, ncomp = 3,
                                   pca.method = c("vertical", "horizontal", "joint"),
                                   lambda = 0, cal.fraction = 0.25,
                                   alpha = 0.1, seed = NULL) {
  if (!inherits(fdataobj, "fdata")) stop("fdataobj must be of class 'fdata'")
  if (!inherits(newdata, "fdata")) stop("newdata must be of class 'fdata'")
  pca.method <- match.arg(pca.method)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  result <- .Call("wrap__conformal_elastic_pcr_rust",
                  fdataobj$data, as.numeric(y), newdata$data,
                  as.numeric(fdataobj$argvals),
                  as.integer(ncomp), pca.method,
                  as.numeric(lambda),
                  as.numeric(cal.fraction), as.numeric(alpha),
                  as.integer(seed))

  if (is.null(result)) stop("conformal.elastic.pcr failed")
  result
}

#' Conformal Classification
#'
#' Constructs conformal prediction sets for functional classification.
#'
#' @param fdataobj An object of class 'fdata' (training data).
#' @param y Class labels (integer vector, 0-indexed).
#' @param newdata An object of class 'fdata' (test data).
#' @param covariates.train Optional scalar covariates for training.
#' @param covariates.test Optional scalar covariates for test.
#' @param ncomp Number of FPC components (default 3).
#' @param classifier Classifier: "lda", "qda", or "knn" (default "lda").
#' @param k.nn k for kNN classifier (default 5).
#' @param score.type Nonconformity score: "lac" or "aps" (default "lac").
#' @param cal.fraction Calibration fraction (default 0.25).
#' @param alpha Miscoverage level (default 0.1).
#' @param seed Random seed.
#'
#' @return A list with components:
#' \describe{
#'   \item{predicted.classes}{Point predictions}
#'   \item{set.sizes}{Size of each prediction set}
#'   \item{average.set.size}{Mean prediction set size}
#'   \item{coverage}{Empirical coverage}
#' }
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- sample(0:1, 50, replace = TRUE)
#' cp <- conformal.classif(fd[1:40, ], y[1:40], fd[41:50, ])
#' }
#'
#' @export
conformal.classif <- function(fdataobj, y, newdata,
                               covariates.train = NULL,
                               covariates.test = NULL,
                               ncomp = 3, classifier = "lda",
                               k.nn = 5, score.type = c("lac", "aps"),
                               cal.fraction = 0.25, alpha = 0.1,
                               seed = NULL) {
  if (!inherits(fdataobj, "fdata")) stop("fdataobj must be of class 'fdata'")
  if (!inherits(newdata, "fdata")) stop("newdata must be of class 'fdata'")
  score.type <- match.arg(score.type)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  ct <- if (!is.null(covariates.train)) as.matrix(covariates.train) else NULL
  cte <- if (!is.null(covariates.test)) as.matrix(covariates.test) else NULL

  result <- .Call("wrap__conformal_classif_rust",
                  fdataobj$data, as.integer(y), newdata$data,
                  ct, cte, as.integer(ncomp), classifier,
                  as.integer(k.nn), score.type,
                  as.numeric(cal.fraction), as.numeric(alpha),
                  as.integer(seed))

  if (is.null(result)) stop("conformal.classif failed")
  result
}

#' Conformal Prediction for Logistic Regression
#'
#' @param fdataobj An object of class 'fdata' (training data).
#' @param y Binary response (0/1).
#' @param newdata An object of class 'fdata' (test data).
#' @param scalar.train Optional scalar covariates for training.
#' @param scalar.test Optional scalar covariates for test.
#' @param ncomp Number of FPC components (default 3).
#' @param max.iter Maximum iterations (default 100).
#' @param tol Convergence tolerance (default 1e-6).
#' @param score.type Nonconformity score: "lac" or "aps".
#' @param cal.fraction Calibration fraction (default 0.25).
#' @param alpha Miscoverage level (default 0.1).
#' @param seed Random seed.
#'
#' @return Same as \code{\link{conformal.classif}}.
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rbinom(50, 1, 0.5)
#' cp <- conformal.logistic(fd[1:40, ], y[1:40], fd[41:50, ])
#' }
#'
#' @export
conformal.logistic <- function(fdataobj, y, newdata,
                                scalar.train = NULL, scalar.test = NULL,
                                ncomp = 3, max.iter = 100, tol = 1e-6,
                                score.type = c("lac", "aps"),
                                cal.fraction = 0.25, alpha = 0.1,
                                seed = NULL) {
  if (!inherits(fdataobj, "fdata")) stop("fdataobj must be of class 'fdata'")
  if (!inherits(newdata, "fdata")) stop("newdata must be of class 'fdata'")
  score.type <- match.arg(score.type)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  st <- if (!is.null(scalar.train)) as.matrix(scalar.train) else NULL
  ste <- if (!is.null(scalar.test)) as.matrix(scalar.test) else NULL

  result <- .Call("wrap__conformal_logistic_rust",
                  fdataobj$data, as.numeric(y), newdata$data,
                  st, ste, as.integer(ncomp),
                  as.integer(max.iter), as.numeric(tol),
                  score.type,
                  as.numeric(cal.fraction), as.numeric(alpha),
                  as.integer(seed))

  if (is.null(result)) stop("conformal.logistic failed")
  result
}

#' Conformal Prediction for Elastic Logistic Regression
#'
#' @param fdataobj An object of class 'fdata' (training data).
#' @param y Binary response (0/1).
#' @param newdata An object of class 'fdata' (test data).
#' @param lambda Regularization parameter (default 0).
#' @param score.type Nonconformity score: "lac" or "aps".
#' @param cal.fraction Calibration fraction (default 0.25).
#' @param alpha Miscoverage level (default 0.1).
#' @param seed Random seed.
#'
#' @return Same as \code{\link{conformal.classif}}.
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- sample(0:1, 50, replace = TRUE)
#' cp <- conformal.elastic.logistic(fd[1:40, ], y[1:40], fd[41:50, ])
#' }
#'
#' @export
conformal.elastic.logistic <- function(fdataobj, y, newdata,
                                        lambda = 0,
                                        score.type = c("lac", "aps"),
                                        cal.fraction = 0.25, alpha = 0.1,
                                        seed = NULL) {
  if (!inherits(fdataobj, "fdata")) stop("fdataobj must be of class 'fdata'")
  if (!inherits(newdata, "fdata")) stop("newdata must be of class 'fdata'")
  score.type <- match.arg(score.type)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  result <- .Call("wrap__conformal_elastic_logistic_rust",
                  fdataobj$data, as.numeric(y), newdata$data,
                  as.numeric(fdataobj$argvals),
                  as.numeric(lambda), score.type,
                  as.numeric(cal.fraction), as.numeric(alpha),
                  as.integer(seed))

  if (is.null(result)) stop("conformal.elastic.logistic failed")
  result
}

# =============================================================================
# Phase C1: CV-Conformal & Jackknife+
# =============================================================================

#' Cross-Conformal (CV+) Regression
#'
#' Constructs prediction intervals using K-fold cross-conformal inference.
#' No data is wasted on a calibration split; each fold provides out-of-fold
#' calibration residuals.
#'
#' @param fdataobj An object of class 'fdata' (training data).
#' @param y Response vector (training).
#' @param newdata An object of class 'fdata' (test data).
#' @param method Regression method: "fregre.lm" or "fregre.np".
#' @param scalar.train Optional scalar covariates for training.
#' @param scalar.test Optional scalar covariates for test.
#' @param ncomp Number of FPC components (default 3, for fregre.lm).
#' @param h.func Functional bandwidth (default 0 = auto, for fregre.np).
#' @param h.scalar Scalar bandwidth (default 0 = auto, for fregre.np).
#' @param n.folds Number of folds (default 5).
#' @param alpha Miscoverage level (default 0.1 for 90 percent intervals).
#' @param seed Random seed.
#'
#' @return A list with components:
#' \describe{
#'   \item{predictions}{Point predictions for test data}
#'   \item{lower}{Lower bounds of prediction intervals}
#'   \item{upper}{Upper bounds of prediction intervals}
#'   \item{residual.quantile}{Calibration residual quantile}
#'   \item{coverage}{Empirical coverage on calibration set}
#' }
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' cp <- cv.conformal.regression(fd[1:40, ], y[1:40], fd[41:50, ])
#' }
#'
#' @export
cv.conformal.regression <- function(fdataobj, y, newdata,
                                     method = c("fregre.lm", "fregre.np"),
                                     scalar.train = NULL, scalar.test = NULL,
                                     ncomp = 3, h.func = 0, h.scalar = 0,
                                     n.folds = 5, alpha = 0.1, seed = NULL) {
  if (!inherits(fdataobj, "fdata")) stop("fdataobj must be of class 'fdata'")
  if (!inherits(newdata, "fdata")) stop("newdata must be of class 'fdata'")
  method <- match.arg(method)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  st <- if (!is.null(scalar.train)) as.matrix(scalar.train) else NULL
  ste <- if (!is.null(scalar.test)) as.matrix(scalar.test) else NULL

  if (method == "fregre.lm") {
    result <- .Call("wrap__cv_conformal_fregre_lm_rust",
                    fdataobj$data, as.numeric(y), newdata$data,
                    st, ste, as.integer(ncomp),
                    as.integer(n.folds), as.numeric(alpha),
                    as.integer(seed))
  } else {
    result <- .Call("wrap__cv_conformal_fregre_np_rust",
                    fdataobj$data, as.numeric(y), newdata$data,
                    as.numeric(fdataobj$argvals), st, ste,
                    as.numeric(h.func), as.numeric(h.scalar),
                    as.integer(n.folds), as.numeric(alpha),
                    as.integer(seed))
  }

  if (is.null(result)) stop("cv.conformal.regression failed")
  result
}

#' Cross-Conformal (CV+) Classification
#'
#' Constructs conformal prediction sets using K-fold cross-conformal inference.
#'
#' @param fdataobj An object of class 'fdata' (training data).
#' @param y Class labels (integer vector, 0-indexed).
#' @param newdata An object of class 'fdata' (test data).
#' @param covariates.train Optional scalar covariates for training.
#' @param covariates.test Optional scalar covariates for test.
#' @param ncomp Number of FPC components (default 3).
#' @param classifier Classifier: "lda", "qda", or "knn" (default "lda").
#' @param k.nn k for kNN classifier (default 5).
#' @param score.type Nonconformity score: "lac" or "aps" (default "lac").
#' @param n.folds Number of folds (default 5).
#' @param alpha Miscoverage level (default 0.1).
#' @param seed Random seed.
#'
#' @return A list with components:
#' \describe{
#'   \item{predicted.classes}{Point predictions}
#'   \item{set.sizes}{Size of each prediction set}
#'   \item{average.set.size}{Mean prediction set size}
#'   \item{coverage}{Empirical coverage}
#' }
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- sample(0:1, 50, replace = TRUE)
#' cp <- cv.conformal.classification(fd[1:40, ], y[1:40], fd[41:50, ])
#' }
#'
#' @export
cv.conformal.classification <- function(fdataobj, y, newdata,
                                          covariates.train = NULL,
                                          covariates.test = NULL,
                                          ncomp = 3, classifier = "lda",
                                          k.nn = 5,
                                          score.type = c("lac", "aps"),
                                          n.folds = 5, alpha = 0.1,
                                          seed = NULL) {
  if (!inherits(fdataobj, "fdata")) stop("fdataobj must be of class 'fdata'")
  if (!inherits(newdata, "fdata")) stop("newdata must be of class 'fdata'")
  score.type <- match.arg(score.type)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  ct <- if (!is.null(covariates.train)) as.matrix(covariates.train) else NULL
  cte <- if (!is.null(covariates.test)) as.matrix(covariates.test) else NULL

  result <- .Call("wrap__cv_conformal_classif_rust",
                  fdataobj$data, as.integer(y), newdata$data,
                  ct, cte, as.integer(ncomp), classifier,
                  as.integer(k.nn), score.type,
                  as.integer(n.folds), as.numeric(alpha),
                  as.integer(seed))

  if (is.null(result)) stop("cv.conformal.classification failed")
  result
}

#' Jackknife+ Regression
#'
#' Constructs prediction intervals using jackknife+ (leave-one-out conformal).
#' Most sample-efficient but most expensive method (n refits).
#'
#' @param fdataobj An object of class 'fdata' (training data).
#' @param y Response vector (training).
#' @param newdata An object of class 'fdata' (test data).
#' @param method Regression method: "fregre.lm" or "fregre.np".
#' @param scalar.train Optional scalar covariates for training.
#' @param scalar.test Optional scalar covariates for test.
#' @param ncomp Number of FPC components (default 3, for fregre.lm).
#' @param h.func Functional bandwidth (default 0 = auto, for fregre.np).
#' @param h.scalar Scalar bandwidth (default 0 = auto, for fregre.np).
#' @param alpha Miscoverage level (default 0.1 for 90 percent intervals).
#' @param seed Random seed (unused, kept for API consistency).
#'
#' @return Same as \code{\link{cv.conformal.regression}}.
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' jp <- jackknife.plus(fd[1:40, ], y[1:40], fd[41:50, ])
#' }
#'
#' @export
jackknife.plus <- function(fdataobj, y, newdata,
                            method = c("fregre.lm", "fregre.np"),
                            scalar.train = NULL, scalar.test = NULL,
                            ncomp = 3, h.func = 0, h.scalar = 0,
                            alpha = 0.1, seed = NULL) {
  if (!inherits(fdataobj, "fdata")) stop("fdataobj must be of class 'fdata'")
  if (!inherits(newdata, "fdata")) stop("newdata must be of class 'fdata'")
  method <- match.arg(method)

  st <- if (!is.null(scalar.train)) as.matrix(scalar.train) else NULL
  ste <- if (!is.null(scalar.test)) as.matrix(scalar.test) else NULL

  if (method == "fregre.lm") {
    result <- .Call("wrap__jackknife_plus_fregre_lm_rust",
                    fdataobj$data, as.numeric(y), newdata$data,
                    st, ste, as.integer(ncomp),
                    as.numeric(alpha))
  } else {
    result <- .Call("wrap__jackknife_plus_fregre_np_rust",
                    fdataobj$data, as.numeric(y), newdata$data,
                    as.numeric(fdataobj$argvals), st, ste,
                    as.numeric(h.func), as.numeric(h.scalar),
                    as.numeric(alpha))
  }

  if (is.null(result)) stop("jackknife.plus failed")
  result
}

# =============================================================================
# Phase C2: Conformal Generic (trait-based)
# =============================================================================

#' Generic Conformal Regression
#'
#' Split-conformal prediction intervals using a pre-fitted fregre.lm model.
#' Uses the model's FPCA components for fast prediction without refitting.
#'
#' @param model A fitted \code{fregre.lm} model object.
#' @param fdataobj An object of class 'fdata' (training data).
#' @param y Response vector (training).
#' @param newdata An object of class 'fdata' (test data).
#' @param scalar.train Optional scalar covariates for training.
#' @param scalar.test Optional scalar covariates for test.
#' @param cal.fraction Fraction of data for calibration (default 0.25).
#' @param alpha Miscoverage level (default 0.1).
#' @param seed Random seed.
#'
#' @return Same as \code{\link{conformal.fregre.lm}}.
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' model <- fregre.lm(fd, y)
#' cp <- conformal.generic.regression(model, fd, y, fd[1:10, ])
#' }
#'
#' @export
conformal.generic.regression <- function(model, fdataobj, y, newdata,
                                           scalar.train = NULL,
                                           scalar.test = NULL,
                                           cal.fraction = 0.25,
                                           alpha = 0.1, seed = NULL) {
  if (!inherits(fdataobj, "fdata")) stop("fdataobj must be of class 'fdata'")
  if (!inherits(newdata, "fdata")) stop("newdata must be of class 'fdata'")
  if (!inherits(model, "fregre.lm")) stop("model must be of class 'fregre.lm'")
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  st <- if (!is.null(scalar.train)) as.matrix(scalar.train) else NULL
  ste <- if (!is.null(scalar.test)) as.matrix(scalar.test) else NULL

  # Extract FPCA components from the model
  fpca_mean <- as.numeric(model$.fpca_mean)
  fpca_rotation <- model$.fpca_rotation
  fpca_scores <- model$.fpca_scores
  ncomp <- model$ncomp
  coefficients <- as.numeric(model$coefficients)
  gamma <- if (!is.null(model$gamma)) as.numeric(model$gamma) else numeric(0)

  result <- .Call("wrap__conformal_generic_regression_lm_rust",
                  fdataobj$data, as.numeric(y), newdata$data,
                  st, ste, fpca_mean, fpca_rotation, fpca_scores,
                  as.integer(ncomp), coefficients, gamma,
                  as.numeric(cal.fraction), as.numeric(alpha),
                  as.integer(seed))

  if (is.null(result)) stop("conformal.generic.regression failed")
  result
}

#' Generic Conformal Classification
#'
#' Split-conformal prediction sets using a pre-fitted functional.logistic model.
#'
#' @param model A fitted model object from \code{\link{functional.logistic}}.
#' @param fdataobj An object of class 'fdata' (training data).
#' @param y Binary response (0/1).
#' @param newdata An object of class 'fdata' (test data).
#' @param scalar.train Optional scalar covariates for training.
#' @param scalar.test Optional scalar covariates for test.
#' @param score.type Nonconformity score: "lac" or "aps" (default "lac").
#' @param cal.fraction Calibration fraction (default 0.25).
#' @param alpha Miscoverage level (default 0.1).
#' @param seed Random seed.
#'
#' @return Same as \code{\link{conformal.classif}}.
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rbinom(50, 1, 0.5)
#' model <- functional.logistic(fd, y)
#' cp <- conformal.generic.classification(model, fd, y, fd[1:10, ])
#' }
#'
#' @export
conformal.generic.classification <- function(model, fdataobj, y, newdata,
                                               scalar.train = NULL,
                                               scalar.test = NULL,
                                               score.type = c("lac", "aps"),
                                               cal.fraction = 0.25,
                                               alpha = 0.1, seed = NULL) {
  if (!inherits(fdataobj, "fdata")) stop("fdataobj must be of class 'fdata'")
  if (!inherits(newdata, "fdata")) stop("newdata must be of class 'fdata'")
  if (!inherits(model, "fregre.logistic")) {
    stop("model must be of class 'fregre.logistic'")
  }
  score.type <- match.arg(score.type)
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)

  st <- if (!is.null(scalar.train)) as.matrix(scalar.train) else NULL
  ste <- if (!is.null(scalar.test)) as.matrix(scalar.test) else NULL

  # Extract FPCA components from the model
  fpca_mean <- as.numeric(model$.fpca_mean)
  fpca_rotation <- model$.fpca_rotation
  fpca_scores <- model$.fpca_scores
  ncomp <- model$ncomp
  intercept <- model$intercept
  coefficients <- as.numeric(model$coefficients)
  gamma <- if (!is.null(model$gamma)) as.numeric(model$gamma) else numeric(0)

  result <- .Call("wrap__conformal_generic_classification_logistic_rust",
                  fdataobj$data, as.integer(y), newdata$data,
                  st, ste, fpca_mean, fpca_rotation, fpca_scores,
                  as.integer(ncomp), as.numeric(intercept),
                  coefficients, gamma,
                  score.type,
                  as.numeric(cal.fraction), as.numeric(alpha),
                  as.integer(seed))

  if (is.null(result)) stop("conformal.generic.classification failed")
  result
}
