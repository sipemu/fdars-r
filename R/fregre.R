#' Functional Regression
#'
#' Functions for functional regression models.

#' Functional Principal Component Regression
#'
#' Fits a functional linear model using principal component regression.
#'
#' @param fdataobj An object of class 'fdata' (functional covariate).
#' @param y Response vector.
#' @param ncomp Number of principal components to use.
#' @param ... Additional arguments.
#'
#' @return A fitted regression object of class 'fregre.fd' with components:
#'   \item{coefficients}{Beta coefficient function values}
#'   \item{intercept}{Intercept term}
#'   \item{fitted.values}{Fitted values}
#'   \item{residuals}{Residuals}
#'   \item{ncomp}{Number of components used}
#'   \item{mean.X}{Mean of functional covariate (for prediction)}
#'   \item{mean.y}{Mean of response (for prediction)}
#'   \item{rotation}{PC loadings (for prediction)}
#'   \item{l}{Indices of selected components}
#'   \item{lm}{Underlying linear model}
#'   \item{sr2}{Residual variance}
#'   \item{fdataobj}{Original functional data}
#'   \item{y}{Response vector}
#'   \item{call}{The function call}
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' result <- fregre.pc(fd, y, ncomp = 3)
#' result$ncomp
#' }
#'
#' @export
fregre.pc <- function(fdataobj, y, ncomp = NULL, ...) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  if (isTRUE(fdataobj$fdata2d)) {
    stop("fregre.pc not yet implemented for 2D functional data")
  }

  n <- nrow(fdataobj$data)
  if (length(y) != n) {
    stop("Length of y must equal number of curves in fdataobj")
  }

  # Store means for prediction
  X_mean <- colMeans(fdataobj$data)
  y_mean <- mean(y)

  # Center the functional data
  X_centered <- scale(fdataobj$data, center = TRUE, scale = FALSE)
  y_centered <- y - y_mean

  # Compute SVD
  svd_result <- svd(X_centered)

  # Determine number of components
  if (is.null(ncomp)) {
    # Use enough components to explain 95% variance
    var_explained <- cumsum(svd_result$d^2) / sum(svd_result$d^2)
    ncomp <- min(which(var_explained >= 0.95), n - 1)
    ncomp <- max(ncomp, 1)
  }

  ncomp <- min(ncomp, length(svd_result$d))
  l <- seq_len(ncomp)

  # Get PC scores (rotation/loadings)
  rotation <- svd_result$v[, l, drop = FALSE]
  scores <- X_centered %*% rotation
  colnames(scores) <- paste0("PC", l)

  # Fit OLS on scores
  scores_df <- as.data.frame(scores)
  lm_fit <- lm(y_centered ~ ., data = scores_df)

  # Compute beta coefficient function
  beta_coef <- rotation %*% coef(lm_fit)[-1]  # Exclude intercept

  # Compute intercept and fitted values
  intercept <- as.numeric(y_mean - sum(X_mean * beta_coef))
  fitted_values <- as.vector(fdataobj$data %*% beta_coef) + intercept

  # Residual variance
  residuals <- y - fitted_values
  sr2 <- sum(residuals^2) / (n - ncomp - 1)

  structure(
    list(
      coefficients = beta_coef,
      intercept = intercept,
      fitted.values = fitted_values,
      residuals = residuals,
      ncomp = ncomp,
      mean.X = X_mean,
      mean.y = y_mean,
      rotation = rotation,
      l = l,
      lm = lm_fit,
      sr2 = sr2,
      svd = svd_result,
      fdataobj = fdataobj,
      y = y,
      call = match.call()
    ),
    class = "fregre.fd"
  )
}

#' Functional Basis Regression
#'
#' Fits a functional linear model using basis expansion (ridge regression).
#' Uses the anofox-regression Rust backend for efficient L2-regularized regression.
#'
#' @param fdataobj An object of class 'fdata' (functional covariate).
#' @param y Response vector.
#' @param basis.x Basis for the functional covariate (currently ignored).
#' @param basis.b Basis for the coefficient function (currently ignored).
#' @param lambda Smoothing/regularization parameter (L2 penalty).
#' @param ... Additional arguments.
#'
#' @return A fitted regression object of class 'fregre.fd' with components:
#'   \item{coefficients}{Beta coefficient function values}
#'   \item{intercept}{Intercept term}
#'   \item{fitted.values}{Fitted values}
#'   \item{residuals}{Residuals}
#'   \item{lambda}{Regularization parameter used}
#'   \item{r.squared}{R-squared (coefficient of determination)}
#'   \item{mean.X}{Mean of functional covariate (for prediction)}
#'   \item{mean.y}{Mean of response (for prediction)}
#'   \item{sr2}{Residual variance}
#'   \item{fdataobj}{Original functional data}
#'   \item{y}{Response vector}
#'   \item{call}{The function call}
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' result <- fregre.basis(fd, y, lambda = 0.1)
#' result$r.squared
#' }
#'
#' @export
fregre.basis <- function(fdataobj, y, basis.x = NULL, basis.b = NULL,
                         lambda = 0, ...) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  n <- nrow(fdataobj$data)
  m <- ncol(fdataobj$data)

  if (length(y) != n) {
    stop("Length of y must equal number of curves in fdataobj")
  }

  # Store means for prediction
  X_mean <- colMeans(fdataobj$data)
  y_mean <- mean(y)

  # Note: Rust ridge regression backend is temporarily disabled due to CRAN

  # compatibility constraints (requires faer 0.23+ which needs Rust 1.84+)
  # Always use the pure R implementation for now
  .fregre_basis_r(fdataobj, y, lambda, X_mean, y_mean)
}

# Internal R implementation for underdetermined systems
.fregre_basis_r <- function(fdataobj, y, lambda, X_mean, y_mean) {
  n <- nrow(fdataobj$data)
  m <- ncol(fdataobj$data)

  # Center the data
  X_centered <- scale(fdataobj$data, center = TRUE, scale = FALSE)
  y_centered <- y - y_mean

  # Ridge regression on centered data: beta = (X'X + lambda*I)^-1 X'y
  XtX <- crossprod(X_centered)
  Xty <- crossprod(X_centered, y_centered)

  if (lambda > 0) {
    XtX <- XtX + lambda * diag(m)
  }

  beta <- solve(XtX, Xty)

  # Compute intercept and fitted values
  intercept <- as.numeric(y_mean - sum(X_mean * beta))
  fitted <- as.vector(fdataobj$data %*% beta) + intercept
  residuals <- y - fitted

  # Compute R-squared
  ss_tot <- sum((y - y_mean)^2)
  ss_res <- sum(residuals^2)
  r_squared <- if (ss_tot > 0) 1 - ss_res / ss_tot else 0

  # Residual variance
  df <- n - m - 1
  if (df <= 0) df <- 1
  sr2 <- sum(residuals^2) / df

  structure(
    list(
      coefficients = beta,
      intercept = intercept,
      fitted.values = fitted,
      residuals = residuals,
      lambda = lambda,
      r.squared = r_squared,
      mean.X = X_mean,
      mean.y = y_mean,
      sr2 = sr2,
      fdataobj = fdataobj,
      y = y,
      call = match.call()
    ),
    class = "fregre.fd"
  )
}

#' Nonparametric Functional Regression
#'
#' Fits a functional regression model using kernel smoothing (Nadaraya-Watson).
#' Supports fixed bandwidth (h), or k-nearest neighbors with global (kNN.gCV)
#' or local (kNN.lCV) cross-validation.
#'
#' @param fdataobj An object of class 'fdata' (functional covariate).
#' @param y Response vector.
#' @param h Bandwidth parameter. If NULL and knn is NULL, computed automatically.
#' @param knn Number of nearest neighbors to consider for bandwidth selection.
#'   Only used when \code{type.S} is "kNN.gCV" or "kNN.lCV".
#' @param type.S Type of smoother: "S.NW" for Nadaraya-Watson with fixed h (default),
#'   "kNN.gCV" for k-NN with global CV (single k for all observations),
#'   "kNN.lCV" for k-NN with local CV (different k per observation).
#' @param Ker Kernel type for smoothing. Default is "norm" (Gaussian).
#' @param metric Distance metric function. Default is metric.lp.
#' @param ... Additional arguments passed to metric function.
#'
#' @return A fitted regression object of class 'fregre.np' with components:
#' \describe{
#'   \item{fitted.values}{Fitted values}
#'   \item{residuals}{Residuals}
#'   \item{h.opt}{Optimal/used bandwidth (for type.S = "S.NW")}
#'   \item{knn}{Number of neighbors used (for kNN methods)}
#'   \item{k.opt}{Optimal k value(s) - scalar for global, vector for local}
#'   \item{type.S}{Type of smoother used}
#'   \item{Ker}{Kernel type used}
#'   \item{fdataobj}{Original functional data}
#'   \item{y}{Response vector}
#'   \item{mdist}{Distance matrix}
#'   \item{sr2}{Residual variance}
#'   \item{metric}{Metric function used}
#'   \item{call}{The function call}
#' }
#'
#' @details
#' Three smoothing approaches are available:
#'
#' \strong{Fixed bandwidth (type.S = "S.NW")}:
#' Uses a single bandwidth h for all predictions. If h is not provided,
#' it is set to the median of non-zero pairwise distances.
#'
#' \strong{k-NN Global CV (type.S = "kNN.gCV")}:
#' Selects a single optimal k for all observations using leave-one-out
#' cross-validation. The bandwidth at each point is set to include k neighbors.
#'
#' \strong{k-NN Local CV (type.S = "kNN.lCV")}:
#' Selects an optimal k_i for each observation i, allowing adaptive smoothing.
#' Useful when the data has varying density across the functional space.
#'
#' @export
#' @examples
#' # Create functional data
#' t <- seq(0, 1, length.out = 50)
#' n <- 50
#' X <- matrix(0, n, 50)
#' for (i in 1:n) X[i, ] <- sin(2*pi*t) * i/n + rnorm(50, sd = 0.1)
#' y <- rowMeans(X) + rnorm(n, sd = 0.1)
#' fd <- fdata(X, argvals = t)
#'
#' # Fixed bandwidth
#' fit1 <- fregre.np(fd, y, h = 0.5)
#'
#' # k-NN with global CV
#' fit2 <- fregre.np(fd, y, type.S = "kNN.gCV", knn = 20)
#'
#' # k-NN with local CV
#' fit3 <- fregre.np(fd, y, type.S = "kNN.lCV", knn = 20)
fregre.np <- function(fdataobj, y, h = NULL, knn = NULL,
                      type.S = c("S.NW", "kNN.gCV", "kNN.lCV"),
                      Ker = "norm", metric = metric.lp, ...) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  type.S <- match.arg(type.S)
  n <- nrow(fdataobj$data)

  if (length(y) != n) {
    stop("Length of y must equal number of curves in fdataobj")
  }

  # Compute distance matrix
  D <- as.matrix(metric(fdataobj, ...))

  result <- list(
    fdataobj = fdataobj,
    y = y,
    mdist = D,
    Ker = Ker,
    type.S = type.S,
    metric = metric,
    call = match.call()
  )

  if (type.S == "S.NW") {
    # Fixed bandwidth Nadaraya-Watson
    if (is.null(h)) {
      d_vec <- D[lower.tri(D)]
      h <- median(d_vec[d_vec > 0])
    }

    H <- matrix(0, n, n)
    fitted <- numeric(n)

    for (i in seq_len(n)) {
      weights <- exp(-0.5 * (D[i, ] / h)^2)
      weights[i] <- 0  # Leave-one-out
      if (sum(weights) > 0) {
        weights <- weights / sum(weights)
      }
      H[i, ] <- weights
      fitted[i] <- sum(weights * y)
    }

    result$h.opt <- h
    result$fitted.values <- fitted
    result$H <- H

  } else if (type.S == "kNN.gCV") {
    # k-NN with Global CV
    if (is.null(knn)) knn <- min(20L, n - 2)
    knn <- as.integer(min(knn, n - 2))

    cv_result <- .Call("wrap__knn_gcv", D, as.numeric(y), knn)

    result$knn <- knn
    result$k.opt <- cv_result$k_opt
    result$fitted.values <- cv_result$yhat
    result$mse <- cv_result$mse

  } else if (type.S == "kNN.lCV") {
    # k-NN with Local CV
    if (is.null(knn)) knn <- min(20L, n - 2)
    knn <- as.integer(min(knn, n - 2))

    cv_result <- .Call("wrap__knn_lcv", D, as.numeric(y), knn)

    result$knn <- knn
    result$k.opt <- cv_result$k_opt
    result$fitted.values <- cv_result$yhat
    result$mse <- cv_result$mse
  }

  result$residuals <- y - result$fitted.values

  # Residual variance
  df <- n - 1
  if (!is.null(result$H)) {
    df <- n - sum(diag(result$H))
  }
  if (df <= 0) df <- 1
  result$sr2 <- sum(result$residuals^2) / df

  structure(result, class = "fregre.np")
}

#' Print method for fregre objects
#' @param x A fregre.fd object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object \code{x}.
#' @export
print.fregre.fd <- function(x, ...) {
  cat("Functional regression model\n")
  cat("  Number of observations:", length(x$y), "\n")
  cat("  R-squared:", 1 - sum(x$residuals^2) / sum((x$y - mean(x$y))^2), "\n")
  invisible(x)
}

#' Print method for fregre.np objects
#' @param x A fregre.np object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object \code{x}.
#' @export
print.fregre.np <- function(x, ...) {
  cat("Nonparametric functional regression model\n")
  cat("  Number of observations:", length(x$y), "\n")
  cat("  Smoother type:", x$type.S, "\n")
  if (!is.null(x$h.opt)) {
    cat("  Bandwidth:", x$h.opt, "\n")
  }
  if (!is.null(x$k.opt)) {
    if (length(x$k.opt) == 1) {
      cat("  Optimal k:", x$k.opt, "\n")
    } else {
      cat("  Optimal k (local):", "min =", min(x$k.opt), ", max =", max(x$k.opt),
          ", median =", median(x$k.opt), "\n")
    }
  }
  r2 <- 1 - sum(x$residuals^2) / sum((x$y - mean(x$y))^2)
  cat("  R-squared:", round(r2, 4), "\n")
  invisible(x)
}

#' Cross-Validation for Functional PC Regression
#'
#' Performs k-fold cross-validation to select the optimal number of
#' principal components for functional PC regression.
#'
#' @param fdataobj An object of class 'fdata' (functional covariate).
#' @param y Response vector.
#' @param kfold Number of folds for cross-validation (default 10).
#' @param ncomp.range Range of number of components to try.
#'   Default is 1 to min(n-1, ncol(data)).
#' @param seed Random seed for fold assignment.
#' @param ... Additional arguments passed to fregre.pc.
#'
#' @return A list with components:
#' \describe{
#'   \item{optimal.ncomp}{Optimal number of components}
#'   \item{cv.errors}{Mean squared prediction error for each ncomp}
#'   \item{cv.se}{Standard error of cv.errors}
#'   \item{model}{Fitted model with optimal ncomp}
#' }
#'
#' @export
#' @examples
#' # Create functional data with a linear relationship
#' set.seed(42)
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 100, 50)
#' for (i in 1:100) X[i, ] <- sin(2*pi*t) * i/100 + rnorm(50, sd = 0.1)
#' beta_true <- cos(2*pi*t)
#' y <- X %*% beta_true + rnorm(100, sd = 0.5)
#' fd <- fdata(X, argvals = t)
#'
#' # Cross-validate to find optimal number of PCs
#' cv_result <- fregre.pc.cv(fd, y, ncomp.range = 1:10)
fregre.pc.cv <- function(fdataobj, y, kfold = 10, ncomp.range = NULL,
                         seed = NULL, ...) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  n <- nrow(fdataobj$data)
  m <- ncol(fdataobj$data)

  if (length(y) != n) {
    stop("Length of y must equal number of curves in fdataobj")
  }

  # Set default ncomp.range
  if (is.null(ncomp.range)) {
    max_comp <- min(n - 1, m)
    ncomp.range <- seq_len(min(max_comp, 15))
  }
  ncomp.range <- ncomp.range[ncomp.range < n & ncomp.range <= m]

  if (length(ncomp.range) == 0) {
    stop("No valid values in ncomp.range")
  }

  # Create fold assignments
  if (!is.null(seed)) set.seed(seed)
  folds <- sample(rep(seq_len(kfold), length.out = n))

  # Initialize storage for CV errors
  cv_errors <- matrix(0, nrow = length(ncomp.range), ncol = kfold)

  for (k in seq_len(kfold)) {
    # Split data
    test_idx <- which(folds == k)
    train_idx <- which(folds != k)

    # Create train/test fdata objects
    fd_train <- fdataobj[train_idx, ]
    fd_test <- fdataobj[test_idx, ]
    y_train <- y[train_idx]
    y_test <- y[test_idx]

    for (j in seq_along(ncomp.range)) {
      ncomp <- ncomp.range[j]

      # Fit model on training data
      tryCatch({
        fit <- fregre.pc(fd_train, y_train, ncomp = ncomp, ...)

        # Predict on test data
        y_pred <- as.vector(fd_test$data %*% fit$coefficients) + fit$intercept

        # Compute MSE
        cv_errors[j, k] <- mean((y_test - y_pred)^2)
      }, error = function(e) {
        cv_errors[j, k] <<- NA
      })
    }
  }

  # Compute mean and SE of CV errors
  mean_cv <- rowMeans(cv_errors, na.rm = TRUE)
  se_cv <- apply(cv_errors, 1, sd, na.rm = TRUE) / sqrt(kfold)
  names(mean_cv) <- names(se_cv) <- ncomp.range

  # Find optimal ncomp (minimum CV error)
  optimal_idx <- which.min(mean_cv)
  optimal_ncomp <- ncomp.range[optimal_idx]

  # Fit final model with optimal ncomp
  final_model <- fregre.pc(fdataobj, y, ncomp = optimal_ncomp, ...)

  list(
    optimal.ncomp = optimal_ncomp,
    cv.errors = mean_cv,
    cv.se = se_cv,
    model = final_model
  )
}

#' Cross-Validation for Functional Basis Regression
#'
#' Performs k-fold cross-validation to select the optimal regularization
#' parameter (lambda) for functional basis regression.
#'
#' @param fdataobj An object of class 'fdata' (functional covariate).
#' @param y Response vector.
#' @param kfold Number of folds for cross-validation (default 10).
#' @param lambda.range Range of lambda values to try.
#'   Default is 10^seq(-4, 4, length.out = 20).
#' @param seed Random seed for fold assignment.
#' @param ... Additional arguments passed to fregre.basis.
#'
#' @return A list with components:
#' \describe{
#'   \item{optimal.lambda}{Optimal regularization parameter}
#'   \item{cv.errors}{Mean squared prediction error for each lambda}
#'   \item{cv.se}{Standard error of cv.errors}
#'   \item{model}{Fitted model with optimal lambda}
#' }
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' cv_result <- fregre.basis.cv(fd, y, kfold = 5,
#'   lambda.range = 10^seq(-2, 2, length.out = 10))
#' cv_result$optimal.lambda
#' }
#'
#' @export
fregre.basis.cv <- function(fdataobj, y, kfold = 10, lambda.range = NULL,
                            seed = NULL, ...) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  n <- nrow(fdataobj$data)

  if (length(y) != n) {
    stop("Length of y must equal number of curves in fdataobj")
  }

  # Set default lambda.range
  if (is.null(lambda.range)) {
    lambda.range <- 10^seq(-4, 4, length.out = 20)
  }

  # Create fold assignments
  if (!is.null(seed)) set.seed(seed)
  folds <- sample(rep(seq_len(kfold), length.out = n))

  # Initialize storage for CV errors
  cv_errors <- matrix(0, nrow = length(lambda.range), ncol = kfold)

  for (k in seq_len(kfold)) {
    # Split data
    test_idx <- which(folds == k)
    train_idx <- which(folds != k)

    # Create train/test fdata objects
    fd_train <- fdataobj[train_idx, ]
    fd_test <- fdataobj[test_idx, ]
    y_train <- y[train_idx]
    y_test <- y[test_idx]

    for (j in seq_along(lambda.range)) {
      lambda <- lambda.range[j]

      # Fit model on training data
      tryCatch({
        fit <- fregre.basis(fd_train, y_train, lambda = lambda, ...)

        # Predict on test data
        y_pred <- as.vector(fd_test$data %*% fit$coefficients)

        # Compute MSE
        cv_errors[j, k] <- mean((y_test - y_pred)^2)
      }, error = function(e) {
        cv_errors[j, k] <<- NA
      })
    }
  }

  # Compute mean and SE of CV errors
  mean_cv <- rowMeans(cv_errors, na.rm = TRUE)
  se_cv <- apply(cv_errors, 1, sd, na.rm = TRUE) / sqrt(kfold)
  names(mean_cv) <- names(se_cv) <- lambda.range

  # Find optimal lambda (minimum CV error)
  optimal_idx <- which.min(mean_cv)
  optimal_lambda <- lambda.range[optimal_idx]

  # Fit final model with optimal lambda
  final_model <- fregre.basis(fdataobj, y, lambda = optimal_lambda, ...)

  list(
    optimal.lambda = optimal_lambda,
    cv.errors = mean_cv,
    cv.se = se_cv,
    model = final_model
  )
}

#' Cross-Validation for Nonparametric Functional Regression
#'
#' Performs k-fold cross-validation to select the optimal bandwidth
#' parameter (h) for nonparametric functional regression.
#'
#' @param fdataobj An object of class 'fdata' (functional covariate).
#' @param y Response vector.
#' @param kfold Number of folds for cross-validation (default 10).
#' @param h.range Range of bandwidth values to try. If NULL, automatically
#'   determined from the distance matrix.
#' @param metric Distance metric function. Default is metric.lp.
#' @param seed Random seed for fold assignment.
#' @param ... Additional arguments passed to the metric function.
#'
#' @return A list with components:
#' \describe{
#'   \item{optimal.h}{Optimal bandwidth parameter}
#'   \item{cv.errors}{Mean squared prediction error for each h}
#'   \item{cv.se}{Standard error of cv.errors}
#'   \item{model}{Fitted model with optimal h}
#' }
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' cv_result <- fregre.np.cv(fd, y, kfold = 5)
#' cv_result$optimal.h
#' }
#'
#' @export
fregre.np.cv <- function(fdataobj, y, kfold = 10, h.range = NULL,
                         metric = metric.lp, seed = NULL, ...) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  n <- nrow(fdataobj$data)

  if (length(y) != n) {
    stop("Length of y must equal number of curves in fdataobj")
  }

  # Compute full distance matrix once
  D <- metric(fdataobj, ...)

  # Set default h.range based on distances
  if (is.null(h.range)) {
    d_vec <- D[lower.tri(D)]
    d_vec <- d_vec[d_vec > 0]
    h.range <- quantile(d_vec, probs = seq(0.05, 0.95, length.out = 20))
  }

  # Create fold assignments
  if (!is.null(seed)) set.seed(seed)
  folds <- sample(rep(seq_len(kfold), length.out = n))

  # Initialize storage for CV errors
  cv_errors <- matrix(0, nrow = length(h.range), ncol = kfold)

  for (k in seq_len(kfold)) {
    # Split data
    test_idx <- which(folds == k)
    train_idx <- which(folds != k)

    y_train <- y[train_idx]
    y_test <- y[test_idx]

    # Extract sub-distance matrices
    D_train <- D[train_idx, train_idx]
    D_test_train <- D[test_idx, train_idx, drop = FALSE]

    for (j in seq_along(h.range)) {
      h <- h.range[j]

      tryCatch({
        # Nadaraya-Watson prediction for test set
        n_test <- length(test_idx)
        y_pred <- numeric(n_test)

        for (i in seq_len(n_test)) {
          # Gaussian kernel weights
          weights <- exp(-0.5 * (D_test_train[i, ] / h)^2)
          weights <- weights / sum(weights)
          y_pred[i] <- sum(weights * y_train)
        }

        # Compute MSE
        cv_errors[j, k] <- mean((y_test - y_pred)^2)
      }, error = function(e) {
        cv_errors[j, k] <<- NA
      })
    }
  }

  # Compute mean and SE of CV errors
  mean_cv <- rowMeans(cv_errors, na.rm = TRUE)
  se_cv <- apply(cv_errors, 1, sd, na.rm = TRUE) / sqrt(kfold)
  names(mean_cv) <- names(se_cv) <- h.range

  # Find optimal h (minimum CV error)
  optimal_idx <- which.min(mean_cv)
  optimal_h <- h.range[optimal_idx]

  # Fit final model with optimal h
  final_model <- fregre.np(fdataobj, y, h = optimal_h, metric = metric, ...)

  list(
    optimal.h = optimal_h,
    cv.errors = mean_cv,
    cv.se = se_cv,
    model = final_model
  )
}

# =============================================================================
# Predict Methods
# =============================================================================

#' Predict Method for Functional Regression (fregre.fd)
#'
#' Predictions from a fitted functional regression model (fregre.pc or fregre.basis).
#'
#' @param object A fitted model object of class 'fregre.fd'.
#' @param newdata An fdata object containing new functional data for prediction.
#'   If NULL, returns fitted values from training data.
#' @param ... Additional arguments (ignored).
#'
#' @return A numeric vector of predicted values.
#'
#' @export
#' @examples
#' # Create functional data
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 30, 50)
#' for (i in 1:30) X[i, ] <- sin(2*pi*t) * i/30 + rnorm(50, sd = 0.1)
#' y <- rowMeans(X) + rnorm(30, sd = 0.1)
#' fd <- fdata(X, argvals = t)
#'
#' # Fit model
#' fit <- fregre.pc(fd, y, ncomp = 3)
#'
#' # Predict on new data
#' X_new <- matrix(sin(2*pi*t) * 0.5, nrow = 1)
#' fd_new <- fdata(X_new, argvals = t)
#' predict(fit, fd_new)
predict.fregre.fd <- function(object, newdata = NULL, ...) {
  if (is.null(object)) {
    stop("No fregre.fd object provided")
  }

  # If no new data, return fitted values
  if (is.null(newdata)) {
    return(object$fitted.values)
  }

  # Convert to fdata if needed
  if (!inherits(newdata, "fdata")) {
    newdata <- fdata(newdata,
                     argvals = object$fdataobj$argvals,
                     rangeval = object$fdataobj$rangeval)
  }

  # Get new data matrix
  X_new <- newdata$data
  nn <- nrow(X_new)

  # Check dimensions match
  if (ncol(X_new) != ncol(object$fdataobj$data)) {
    stop("Number of evaluation points in newdata must match training data")
  }

  # Compute predictions using stored coefficients and intercept
  if (!is.null(object$intercept)) {
    y_pred <- as.vector(X_new %*% object$coefficients) + object$intercept
  } else {
    y_pred <- as.vector(X_new %*% object$coefficients)
  }

  names(y_pred) <- rownames(X_new)
  y_pred
}

#' Predict Method for Nonparametric Functional Regression (fregre.np)
#'
#' Predictions from a fitted nonparametric functional regression model.
#'
#' @param object A fitted model object of class 'fregre.np'.
#' @param newdata An fdata object containing new functional data for prediction.
#'   If NULL, returns fitted values from training data.
#' @param ... Additional arguments (ignored).
#'
#' @return A numeric vector of predicted values.
#'
#' @export
#' @examples
#' # Create functional data
#' t <- seq(0, 1, length.out = 50)
#' X <- matrix(0, 30, 50)
#' for (i in 1:30) X[i, ] <- sin(2*pi*t) * i/30 + rnorm(50, sd = 0.1)
#' y <- rowMeans(X) + rnorm(30, sd = 0.1)
#' fd <- fdata(X, argvals = t)
#'
#' # Fit model
#' fit <- fregre.np(fd, y)
#'
#' # Predict on new data
#' X_new <- matrix(sin(2*pi*t) * 0.5, nrow = 1)
#' fd_new <- fdata(X_new, argvals = t)
#' predict(fit, fd_new)
predict.fregre.np <- function(object, newdata = NULL, ...) {
  if (is.null(object)) {
    stop("No fregre.np object provided")
  }

  # If no new data, return fitted values
  if (is.null(newdata)) {
    return(object$fitted.values)
  }

  # Convert to fdata if needed
  if (!inherits(newdata, "fdata")) {
    newdata <- fdata(newdata,
                     argvals = object$fdataobj$argvals,
                     rangeval = object$fdataobj$rangeval)
  }

  # Check dimensions match
  if (ncol(newdata$data) != ncol(object$fdataobj$data)) {
    stop("Number of evaluation points in newdata must match training data")
  }

  nn <- nrow(newdata$data)
  n_train <- nrow(object$fdataobj$data)
  y_train <- object$y

  # Compute distances from new data to training data
  metric_fn <- object$metric
  D_new <- as.matrix(metric_fn(object$fdataobj, newdata))

  type.S <- if (!is.null(object$type.S)) object$type.S else "S.NW"

  if (type.S == "S.NW") {
    # Fixed bandwidth Nadaraya-Watson prediction
    h <- object$h.opt
    y_pred <- numeric(nn)

    for (i in seq_len(nn)) {
      weights <- exp(-0.5 * (D_new[, i] / h)^2)
      if (sum(weights) > 0) {
        weights <- weights / sum(weights)
      } else {
        weights <- rep(1/n_train, n_train)
      }
      y_pred[i] <- sum(weights * y_train)
    }

  } else if (type.S %in% c("kNN.gCV", "kNN.lCV")) {
    # k-NN prediction using Rust backend
    k_opt <- object$k.opt

    if (length(k_opt) == 1) {
      # Global k - pass as single integer
      y_pred <- .Call("wrap__knn_predict", D_new, as.numeric(y_train),
                      as.integer(k_opt), NULL)
    } else {
      # Local k - pass as vector
      y_pred <- .Call("wrap__knn_predict", D_new, as.numeric(y_train),
                      0L, as.integer(k_opt))
    }
  }

  names(y_pred) <- rownames(newdata$data)
  y_pred
}

#' Nonparametric Regression with Multiple Functional Predictors
#'
#' Fits a nonparametric regression model with multiple functional predictors
#' using a weighted combination of distance matrices.
#'
#' @param fdataobj.list A list of fdata objects (functional predictors).
#'   All must have the same number of observations.
#' @param y Response vector (scalar).
#' @param weights Weights for combining distances. Can be:
#'   \itemize{
#'     \item NULL: Equal weights (1/p for each predictor)
#'     \item Numeric vector: Fixed weights (will be normalized to sum to 1)
#'     \item "cv": Cross-validate to find optimal weights
#'   }
#' @param h Bandwidth for Nadaraya-Watson kernel (optional).
#' @param knn Maximum k for k-NN methods.
#' @param type.S Smoother type: "S.NW", "kNN.gCV", or "kNN.lCV".
#' @param Ker Kernel type (default "norm" for Gaussian).
#' @param metric Distance metric function (default metric.lp).
#' @param cv.grid Grid of weight values for CV (only used if weights = "cv").
#'   Default is seq(0, 1, by = 0.1) for 2 predictors.
#' @param cv.folds Number of folds for weight CV (default 5).
#' @param ... Additional arguments passed to metric function.
#'
#' @return An object of class 'fregre.np.multi' containing:
#' \describe{
#'   \item{fdataobj.list}{List of functional predictors}
#'   \item{y}{Response vector}
#'   \item{weights}{Weights used (or optimized)}
#'   \item{weights.cv}{CV results if weights = "cv"}
#'   \item{fitted.values}{Fitted values}
#'   \item{residuals}{Residuals}
#'   \item{D.list}{List of distance matrices}
#'   \item{D.combined}{Combined distance matrix}
#' }
#'
#' @export
#' @examples
#' # Create two functional predictors
#' set.seed(42)
#' n <- 50
#' m <- 30
#' t_grid <- seq(0, 1, length.out = m)
#'
#' X1 <- matrix(0, n, m)
#' X2 <- matrix(0, n, m)
#' for (i in 1:n) {
#'   X1[i, ] <- sin(2 * pi * t_grid) * i/n + rnorm(m, sd = 0.1)
#'   X2[i, ] <- cos(2 * pi * t_grid) * i/n + rnorm(m, sd = 0.1)
#' }
#' y <- rowMeans(X1) + 0.5 * rowMeans(X2) + rnorm(n, sd = 0.1)
#'
#' fd1 <- fdata(X1, argvals = t_grid)
#' fd2 <- fdata(X2, argvals = t_grid)
#'
#' # Fit with equal weights
#' fit1 <- fregre.np.multi(list(fd1, fd2), y)
#'
#' # Fit with cross-validated weights
#' fit2 <- fregre.np.multi(list(fd1, fd2), y, weights = "cv")
fregre.np.multi <- function(fdataobj.list, y, weights = NULL,
                            h = NULL, knn = NULL,
                            type.S = c("S.NW", "kNN.gCV", "kNN.lCV"),
                            Ker = "norm", metric = metric.lp,
                            cv.grid = NULL, cv.folds = 5, ...) {

  # Validate inputs
  if (!is.list(fdataobj.list)) {
    stop("fdataobj.list must be a list of fdata objects")
  }

  p <- length(fdataobj.list)
  if (p < 1) {
    stop("Need at least one functional predictor")
  }

  # Check all are fdata with same n
  n_vec <- sapply(fdataobj.list, function(fd) {
    if (!inherits(fd, "fdata")) {
      stop("All elements of fdataobj.list must be fdata objects")
    }
    nrow(fd$data)
  })

  if (length(unique(n_vec)) > 1) {
    stop("All functional predictors must have the same number of observations")
  }

  n <- n_vec[1]
  type.S <- match.arg(type.S)

  if (length(y) != n) {
    stop("Length of y must equal number of observations (", n, ")")
  }

  # Compute distance matrices for each predictor
  D.list <- lapply(fdataobj.list, function(fd) {
    as.matrix(metric(fd, ...))
  })

  # Normalize each distance matrix to [0, 1] range for fair comparison
  D.list.norm <- lapply(D.list, function(D) {
    max_d <- max(D)
    if (max_d > 0) D / max_d else D
  })

  # Handle weights
  cv_result <- NULL

  if (is.null(weights)) {
    # Equal weights
    weights <- rep(1 / p, p)

  } else if (is.character(weights) && weights == "cv") {
    # Cross-validate weights
    cv_result <- .cv_weights_multi(D.list.norm, y, type.S, h, knn, Ker,
                                   cv.grid, cv.folds, p)
    weights <- cv_result$optimal.weights

  } else if (is.numeric(weights)) {
    if (length(weights) != p) {
      stop("Length of weights must equal number of predictors (", p, ")")
    }
    # Normalize to sum to 1
    weights <- weights / sum(weights)
  }

  # Combine distance matrices
  D.combined <- matrix(0, n, n)
  for (i in seq_len(p)) {
    D.combined <- D.combined + weights[i] * D.list.norm[[i]]
  }

  # Fit using combined distance
  result <- .fit_np_with_distance(D.combined, y, h, knn, type.S, Ker)

  # Build result object
  result$fdataobj.list <- fdataobj.list
  result$y <- y
  result$weights <- weights
  result$weights.cv <- cv_result
  result$D.list <- D.list
  result$D.combined <- D.combined
  result$metric <- metric
  result$type.S <- type.S
  result$Ker <- Ker
  result$call <- match.call()

  class(result) <- "fregre.np.multi"
  result
}

# Internal: Cross-validate weights for multiple predictors
.cv_weights_multi <- function(D.list.norm, y, type.S, h, knn, Ker,
                              cv.grid, cv.folds, p) {
  n <- length(y)

  # Generate weight grid
  if (is.null(cv.grid)) {
    if (p == 2) {
      # For 2 predictors, simple 1D grid
      w1_grid <- seq(0, 1, by = 0.1)
      weight_grid <- lapply(w1_grid, function(w1) c(w1, 1 - w1))
    } else {
      # For >2 predictors, use deterministic simplex grid
      # Generate evenly spaced points using Halton-like sequence
      n_grid <- 20
      primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29)[seq_len(min(p, 10))]
      weight_grid <- lapply(seq_len(n_grid), function(i) {
        w <- numeric(p)
        for (j in seq_len(p)) {
          # Van der Corput sequence for dimension j
          base <- primes[((j - 1) %% length(primes)) + 1]
          val <- 0; f <- 1; ii <- i
          while (ii > 0) { f <- f / base; val <- val + (ii %% base) * f; ii <- ii %/% base }
          w[j] <- val + 0.01
        }
        w / sum(w)
      })
      # Add equal weights and extreme weights
      weight_grid <- c(weight_grid, list(rep(1/p, p)))
      for (i in seq_len(p)) {
        extreme <- rep(0, p)
        extreme[i] <- 1
        weight_grid <- c(weight_grid, list(extreme))
      }
    }
  } else {
    weight_grid <- cv.grid
  }

  # Create folds
  fold_ids <- sample(rep(seq_len(cv.folds), length.out = n))

  # Evaluate each weight combination
  cv_errors <- sapply(weight_grid, function(w) {
    # Combine distances
    D.combined <- matrix(0, n, n)
    for (i in seq_len(p)) {
      D.combined <- D.combined + w[i] * D.list.norm[[i]]
    }

    # K-fold CV
    fold_errors <- sapply(seq_len(cv.folds), function(fold) {
      test_idx <- which(fold_ids == fold)
      train_idx <- which(fold_ids != fold)

      # Fit on training, predict on test
      D_train <- D.combined[train_idx, train_idx]
      y_train <- y[train_idx]
      D_test <- D.combined[test_idx, train_idx]

      # Simple NW prediction for CV
      if (is.null(h)) {
        d_vec <- D_train[lower.tri(D_train)]
        h_cv <- median(d_vec[d_vec > 0])
      } else {
        h_cv <- h
      }

      y_pred <- sapply(seq_len(nrow(D_test)), function(i) {
        weights_nw <- exp(-0.5 * (D_test[i, ] / h_cv)^2)
        if (sum(weights_nw) > 0) {
          sum(weights_nw * y_train) / sum(weights_nw)
        } else {
          mean(y_train)
        }
      })

      mean((y[test_idx] - y_pred)^2)
    })

    mean(fold_errors)
  })

  # Find best weights
  best_idx <- which.min(cv_errors)
  optimal_weights <- weight_grid[[best_idx]]

  list(
    optimal.weights = optimal_weights,
    cv.errors = cv_errors,
    weight.grid = weight_grid,
    best.idx = best_idx
  )
}

# Internal: Fit NP regression with pre-computed distance matrix
.fit_np_with_distance <- function(D, y, h, knn, type.S, Ker) {
  n <- length(y)

  result <- list()

  if (type.S == "S.NW") {
    # Fixed bandwidth Nadaraya-Watson
    if (is.null(h)) {
      d_vec <- D[lower.tri(D)]
      h <- median(d_vec[d_vec > 0])
    }

    H <- matrix(0, n, n)
    fitted <- numeric(n)

    for (i in seq_len(n)) {
      weights <- exp(-0.5 * (D[i, ] / h)^2)
      weights[i] <- 0  # Leave-one-out
      if (sum(weights) > 0) {
        weights <- weights / sum(weights)
      }
      H[i, ] <- weights
      fitted[i] <- sum(weights * y)
    }

    result$h.opt <- h
    result$fitted.values <- fitted
    result$H <- H

  } else if (type.S == "kNN.gCV") {
    # k-NN with Global CV
    if (is.null(knn)) knn <- min(20L, n - 2)
    knn <- as.integer(min(knn, n - 2))

    cv_result <- .Call("wrap__knn_gcv", D, as.numeric(y), knn)

    result$knn <- knn
    result$k.opt <- cv_result$k_opt
    result$fitted.values <- cv_result$yhat
    result$mse <- cv_result$mse

  } else if (type.S == "kNN.lCV") {
    # k-NN with Local CV
    if (is.null(knn)) knn <- min(20L, n - 2)
    knn <- as.integer(min(knn, n - 2))

    cv_result <- .Call("wrap__knn_lcv", D, as.numeric(y), knn)

    result$knn <- knn
    result$k.opt <- cv_result$k_opt
    result$fitted.values <- cv_result$yhat
    result$mse <- cv_result$mse
  }

  result$residuals <- y - result$fitted.values

  # Residual variance
  df <- n - 1
  if (!is.null(result$H)) {
    df <- n - sum(diag(result$H))
  }
  if (df <= 0) df <- 1
  result$sr2 <- sum(result$residuals^2) / df

  result
}

#' Print method for fregre.np.multi
#' @param x A fregre.np.multi object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object \code{x}.
#' @export
print.fregre.np.multi <- function(x, ...) {
  cat("Nonparametric functional regression with multiple predictors\n")
  cat("=============================================================\n")
  cat("Number of observations:", length(x$y), "\n")
  cat("Number of functional predictors:", length(x$fdataobj.list), "\n")
  cat("Smoother type:", x$type.S, "\n")
  cat("Weights:", paste(round(x$weights, 3), collapse = ", "), "\n")

  if (!is.null(x$weights.cv)) {
    cat("  (cross-validated)\n")
  }

  if (!is.null(x$h.opt)) {
    cat("Bandwidth:", round(x$h.opt, 4), "\n")
  }
  if (!is.null(x$k.opt)) {
    if (length(x$k.opt) == 1) {
      cat("Optimal k:", x$k.opt, "\n")
    } else {
      cat("Optimal k (local): min =", min(x$k.opt), ", max =", max(x$k.opt), "\n")
    }
  }

  r2 <- 1 - sum(x$residuals^2) / sum((x$y - mean(x$y))^2)
  cat("R-squared:", round(r2, 4), "\n")

  invisible(x)
}

#' Predict method for fregre.np.multi
#'
#' @param object Fitted fregre.np.multi object.
#' @param newdata.list List of new fdata objects (same length as original).
#' @param ... Additional arguments.
#'
#' @return Predicted values.
#' @export
predict.fregre.np.multi <- function(object, newdata.list = NULL, ...) {
  if (is.null(newdata.list)) {
    return(object$fitted.values)
  }

  p <- length(object$fdataobj.list)
  if (length(newdata.list) != p) {
    stop("newdata.list must have same length as original fdataobj.list (", p, ")")
  }

  # Compute distances from new data to training data
  n_new <- nrow(newdata.list[[1]]$data)
  n_train <- nrow(object$fdataobj.list[[1]]$data)

  D_new_list <- lapply(seq_len(p), function(i) {
    # Combine new and training data
    combined <- fdata(
      rbind(newdata.list[[i]]$data, object$fdataobj.list[[i]]$data),
      argvals = object$fdataobj.list[[i]]$argvals
    )
    D_full <- as.matrix(object$metric(combined, ...))

    # Extract new-to-train distances and normalize
    D_sub <- D_full[seq_len(n_new), (n_new + 1):(n_new + n_train), drop = FALSE]

    # Normalize using training max
    max_train <- max(object$D.list[[i]])
    if (max_train > 0) D_sub / max_train else D_sub
  })

  # Combine with weights
  D_new <- matrix(0, n_new, n_train)
  for (i in seq_len(p)) {
    D_new <- D_new + object$weights[i] * D_new_list[[i]]
  }

  y_train <- object$y
  type.S <- object$type.S

  if (type.S == "S.NW") {
    h <- object$h.opt
    y_pred <- sapply(seq_len(n_new), function(i) {
      weights <- exp(-0.5 * (D_new[i, ] / h)^2)
      if (sum(weights) > 0) {
        sum(weights * y_train) / sum(weights)
      } else {
        mean(y_train)
      }
    })

  } else if (type.S %in% c("kNN.gCV", "kNN.lCV")) {
    k_opt <- object$k.opt

    if (length(k_opt) == 1) {
      y_pred <- .Call("wrap__knn_predict", D_new, as.numeric(y_train),
                      as.integer(k_opt), NULL)
    } else {
      y_pred <- .Call("wrap__knn_predict", D_new, as.numeric(y_train),
                      0L, as.integer(k_opt))
    }
  }

  y_pred
}

# =============================================================================
# Scalar-on-function regression (Rust-backed via fdars-core)
# =============================================================================

#' Functional Linear Model (FPC-based)
#'
#' Fits a functional linear model using FPC regression with the Rust backend.
#' Supports optional scalar covariates.
#'
#' @param fdataobj An object of class 'fdata' (functional covariate).
#' @param y Response vector (scalar).
#' @param scalar.covariates Optional matrix of scalar covariates (n x p).
#' @param ncomp Number of FPC components to use. If NULL, uses min(n-1, 15).
#'
#' @return A fitted regression object of class 'fregre.lm' with components:
#'   \item{intercept}{Intercept term}
#'   \item{beta.t}{Functional coefficient beta(t) as fdata}
#'   \item{gamma}{Scalar covariate coefficients}
#'   \item{fitted.values}{Fitted values}
#'   \item{residuals}{Residuals}
#'   \item{r.squared}{R-squared}
#'   \item{r.squared.adj}{Adjusted R-squared}
#'   \item{ncomp}{Number of FPC components used}
#'   \item{gcv}{GCV criterion value}
#'
#' @seealso \code{\link{fregre.pc}} for the R-native FPC regression,
#'   \code{\link{fregre.lm.cv}} for cross-validated component selection
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' result <- fregre.lm(fd, y, ncomp = 3)
#' result$r.squared
#' }
#'
#' @export
fregre.lm <- function(fdataobj, y, scalar.covariates = NULL, ncomp = NULL) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  n <- nrow(fdataobj$data)
  m <- ncol(fdataobj$data)

  if (length(y) != n) {
    stop("Length of y must equal number of curves in fdataobj")
  }

  if (is.null(ncomp)) {
    ncomp <- min(n - 1, 15)
  }

  sc <- if (!is.null(scalar.covariates)) as.matrix(scalar.covariates) else NULL

  result <- .Call("wrap__fregre_lm_rust", fdataobj$data, as.numeric(y),
                  sc, as.integer(ncomp))

  if (is.null(result)) {
    stop("fregre.lm failed: check data dimensions and ncomp")
  }

  beta_fd <- fdata(matrix(result$beta_t, nrow = 1), argvals = fdataobj$argvals,
                   names = list(main = "beta(t)", xlab = fdataobj$names$xlab,
                                ylab = "beta"))

  structure(
    list(
      intercept = result$intercept,
      beta.t = beta_fd,
      beta.se = result$beta_se,
      gamma = result$gamma,
      coefficients = result$coefficients,
      fitted.values = result$fitted_values,
      residuals = result$residuals,
      r.squared = result$r_squared,
      r.squared.adj = result$r_squared_adj,
      std.errors = result$std_errors,
      ncomp = result$ncomp,
      residual.se = result$residual_se,
      gcv = result$gcv,
      aic = result$aic,
      bic = result$bic,
      fdataobj = fdataobj,
      y = y,
      scalar.covariates = scalar.covariates,
      .fpca_mean = result$fpca_mean,
      .fpca_rotation = matrix(result$fpca_rotation_data,
                              nrow = result$fpca_rotation_nrow,
                              ncol = result$fpca_rotation_ncol),
      .fpca_scores = if (!is.null(result$fpca_scores_data) && length(result$fpca_scores_data) > 0) {
        matrix(result$fpca_scores_data,
               nrow = result$fpca_scores_nrow,
               ncol = result$fpca_scores_ncol)
      } else {
        matrix(0, nrow = 0, ncol = 0)
      },
      call = match.call()
    ),
    class = "fregre.lm"
  )
}

#' Nonparametric Functional Regression with Mixed Predictors
#'
#' Fits a nonparametric kernel regression with functional and scalar predictors
#' using the Rust backend.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param y Response vector.
#' @param scalar.covariates Optional matrix of scalar covariates.
#' @param h.func Bandwidth for functional distance kernel. If NULL, auto-selected.
#' @param h.scalar Bandwidth for scalar covariates kernel. If NULL, auto-selected.
#'
#' @return A fitted regression object of class 'fregre.np'.
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' scalars <- matrix(rnorm(100), nrow = 50, ncol = 2)
#' result <- fregre.np.mixed(fd, y, scalar.covariates = scalars)
#' result$r.squared
#' }
#'
#' @export
fregre.np.mixed <- function(fdataobj, y, scalar.covariates = NULL,
                             h.func = NULL, h.scalar = NULL) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  n <- nrow(fdataobj$data)
  if (length(y) != n) {
    stop("Length of y must equal number of curves in fdataobj")
  }

  if (is.null(h.func)) h.func <- 0.0
  if (is.null(h.scalar)) h.scalar <- 0.0

  sc <- if (!is.null(scalar.covariates)) as.matrix(scalar.covariates) else NULL

  result <- .Call("wrap__fregre_np_mixed_rust", fdataobj$data, as.numeric(y),
                  as.numeric(fdataobj$argvals), sc,
                  as.numeric(h.func), as.numeric(h.scalar))

  if (is.null(result)) {
    stop("fregre.np.mixed failed: check data dimensions")
  }

  structure(
    list(
      fitted.values = result$fitted_values,
      residuals = result$residuals,
      r.squared = result$r_squared,
      h.func = result$h_func,
      h.scalar = result$h_scalar,
      cv.error = result$cv_error,
      fdataobj = fdataobj,
      y = y,
      scalar.covariates = scalar.covariates,
      call = match.call()
    ),
    class = "fregre.np"
  )
}

#' Functional Logistic Regression
#'
#' Binary logistic regression with functional and optional scalar predictors
#' using FPC projection and IRLS.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param y Binary response vector (0/1).
#' @param scalar.covariates Optional matrix of scalar covariates.
#' @param ncomp Number of FPC components (default 3).
#' @param max.iter Maximum IRLS iterations (default 100).
#' @param tol Convergence tolerance (default 1e-6).
#'
#' @return A fitted object of class 'fregre.logistic' with components:
#'   \item{probabilities}{Predicted probabilities P(Y=1)}
#'   \item{predicted.classes}{Predicted class labels (0 or 1)}
#'   \item{accuracy}{Classification accuracy on training data}
#'   \item{log.likelihood}{Log-likelihood at convergence}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' t_grid <- seq(0, 1, length.out = 30)
#' X <- matrix(0, 40, 30)
#' for (i in 1:40) X[i, ] <- sin(2*pi*t_grid) * (2*(i > 20) - 1) + rnorm(30, sd = 0.2)
#' fd <- fdata(X, argvals = t_grid)
#' y <- as.numeric(1:40 > 20)
#' result <- functional.logistic(fd, y, ncomp = 3)
#' result$accuracy
#' }
#'
#' @export
functional.logistic <- function(fdataobj, y, scalar.covariates = NULL,
                                 ncomp = 3, max.iter = 100, tol = 1e-6) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  n <- nrow(fdataobj$data)
  if (length(y) != n) {
    stop("Length of y must equal number of curves in fdataobj")
  }

  sc <- if (!is.null(scalar.covariates)) as.matrix(scalar.covariates) else NULL

  result <- .Call("wrap__functional_logistic_rust", fdataobj$data, as.numeric(y),
                  sc, as.integer(ncomp), as.integer(max.iter), as.numeric(tol))

  if (is.null(result)) {
    stop("functional.logistic failed: check data dimensions and response")
  }

  beta_fd <- fdata(matrix(result$beta_t, nrow = 1), argvals = fdataobj$argvals,
                   names = list(main = "beta(t)", xlab = fdataobj$names$xlab,
                                ylab = "beta"))

  structure(
    list(
      intercept = result$intercept,
      beta.t = beta_fd,
      beta.se = result$beta_se,
      gamma = result$gamma,
      probabilities = result$probabilities,
      predicted.classes = result$predicted_classes,
      ncomp = result$ncomp,
      accuracy = result$accuracy,
      std.errors = result$std_errors,
      coefficients = result$coefficients,
      log.likelihood = result$log_likelihood,
      aic = result$aic,
      bic = result$bic,
      iterations = result$iterations,
      fdataobj = fdataobj,
      y = y,
      scalar.covariates = scalar.covariates,
      .fpca_mean = result$fpca_mean,
      .fpca_rotation = matrix(result$fpca_rotation_data,
                              nrow = result$fpca_rotation_nrow,
                              ncol = result$fpca_rotation_ncol),
      .fpca_scores = matrix(result$fpca_scores_data,
                            nrow = result$fpca_scores_nrow,
                            ncol = result$fpca_scores_ncol),
      call = match.call()
    ),
    class = "fregre.logistic"
  )
}

#' Predict from Functional Logistic Model
#'
#' @param object A fitted object of class 'fregre.logistic'.
#' @param newdata New functional data (fdata object). If NULL, returns
#'   fitted probabilities.
#' @param new.scalar Optional matrix of new scalar covariates.
#' @param type "response" for probabilities (default) or "class" for classes.
#' @param ... Additional arguments (ignored).
#'
#' @return Predicted probabilities or class labels.
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rbinom(50, 1, 0.5)
#' fit <- functional.logistic(fd, y, ncomp = 3)
#' predict(fit)
#' }
#'
#' @export
predict.fregre.logistic <- function(object, newdata = NULL, new.scalar = NULL,
                                     type = c("response", "class"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    probs <- object$probabilities
  } else {
    if (!inherits(newdata, "fdata")) {
      newdata <- fdata(newdata, argvals = object$fdataobj$argvals)
    }
    X_centered <- sweep(newdata$data, 2, object$.fpca_mean)
    scores <- X_centered %*% object$.fpca_rotation
    ncomp <- object$ncomp
    coefs <- object$coefficients
    fpc_coefs <- coefs[seq_len(ncomp) + 1L]
    eta <- rep(object$intercept, nrow(newdata$data))
    eta <- eta + as.vector(scores[, seq_len(ncomp), drop = FALSE] %*% fpc_coefs)
    if (!is.null(new.scalar) && length(object$gamma) > 0) {
      new.scalar <- as.matrix(new.scalar)
      eta <- eta + as.vector(new.scalar %*% object$gamma)
    }
    probs <- 1 / (1 + exp(-eta))
  }
  if (type == "class") {
    return(as.integer(probs >= 0.5))
  }
  probs
}

#' Model Selection for Number of FPC Components
#'
#' Selects the optimal number of FPC components for scalar-on-function
#' regression using AIC, BIC, or GCV criterion.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param y Response vector (scalar).
#' @param scalar.covariates Optional matrix of scalar covariates.
#' @param max.ncomp Maximum number of components to evaluate (default 15).
#' @param criterion Selection criterion: "aic", "bic", or "gcv".
#'
#' @return A list with components:
#'   \item{best.ncomp}{Optimal number of components}
#'   \item{ncomp}{Vector of tested component counts}
#'   \item{aic}{AIC values}
#'   \item{bic}{BIC values}
#'   \item{gcv}{GCV values}
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' model.selection.ncomp(fd, y, criterion = "bic")
#' }
#'
#' @export
model.selection.ncomp <- function(fdataobj, y, scalar.covariates = NULL,
                                   max.ncomp = 15, criterion = c("aic", "bic", "gcv")) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }
  criterion <- match.arg(criterion)

  sc <- if (!is.null(scalar.covariates)) as.matrix(scalar.covariates) else NULL

  result <- .Call("wrap__model_selection_ncomp_rust", fdataobj$data,
                  as.numeric(y), sc, as.integer(max.ncomp), criterion)

  if (is.null(result)) {
    stop("model.selection.ncomp failed")
  }

  list(
    best.ncomp = result$best_ncomp,
    ncomp = result$ncomp,
    aic = result$aic,
    bic = result$bic,
    gcv = result$gcv
  )
}

#' Cross-Validation for FPC Component Selection (fregre.lm)
#'
#' Uses the Rust backend to select the optimal number of FPC components
#' via k-fold cross-validation for the functional linear model.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param y Response vector.
#' @param scalar.covariates Optional matrix of scalar covariates.
#' @param k.range Range of FPC component counts to try.
#' @param nfold Number of CV folds (default 10).
#'
#' @return A list with \code{optimal.k}, \code{cv.errors}, and \code{model}.
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' cv_result <- fregre.lm.cv(fd, y, k.range = 1:5, nfold = 5)
#' cv_result$optimal.k
#' }
#'
#' @export
fregre.lm.cv <- function(fdataobj, y, scalar.covariates = NULL,
                          k.range = NULL, nfold = 10) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  n <- nrow(fdataobj$data)
  m <- ncol(fdataobj$data)

  if (is.null(k.range)) {
    k.range <- c(1L, min(n - 1L, m, 15L))
  } else {
    k.range <- c(min(k.range), max(k.range))
  }

  sc <- if (!is.null(scalar.covariates)) as.matrix(scalar.covariates) else NULL

  result <- .Call("wrap__fregre_cv_rust", fdataobj$data, as.numeric(y),
                  sc, as.integer(k.range[1]), as.integer(k.range[2]),
                  as.integer(nfold))

  if (is.null(result)) {
    stop("fregre.lm.cv failed")
  }

  final_model <- fregre.lm(fdataobj, y, scalar.covariates, ncomp = result$optimal_k)

  list(
    optimal.k = result$optimal_k,
    cv.errors = result$cv_errors,
    k.values = result$k_values,
    min.cv.error = result$min_cv_error,
    model = final_model
  )
}

#' Print method for fregre.lm objects
#' @param x A fregre.lm object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object \code{x}.
#' @export
print.fregre.lm <- function(x, ...) {
  cat("Functional Linear Model (FPC-based)\n")
  cat("====================================\n")
  cat("  Number of observations:", length(x$y), "\n")
  cat("  Number of FPC components:", x$ncomp, "\n")
  if (length(x$gamma) > 0) {
    cat("  Scalar covariates:", length(x$gamma), "\n")
  }
  cat("  R-squared:", round(x$r.squared, 4), "\n")
  cat("  Adjusted R-squared:", round(x$r.squared.adj, 4), "\n")
  cat("  GCV:", round(x$gcv, 4), "\n")
  invisible(x)
}

#' Print method for fregre.logistic objects
#' @param x A fregre.logistic object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object \code{x}.
#' @export
print.fregre.logistic <- function(x, ...) {
  cat("Functional Logistic Regression\n")
  cat("==============================\n")
  cat("  Number of observations:", length(x$y), "\n")
  cat("  Number of FPC components:", x$ncomp, "\n")
  cat("  Accuracy:", round(x$accuracy, 4), "\n")
  cat("  Log-likelihood:", round(x$log.likelihood, 2), "\n")
  cat("  Iterations:", x$iterations, "\n")
  invisible(x)
}

#' Predict method for fregre.lm objects
#' @param object A fregre.lm object.
#' @param newdata An fdata object (or NULL for fitted values).
#' @param new.scalar Optional new scalar covariates matrix.
#' @param ... Additional arguments (ignored).
#' @return Predicted values.
#' @export
predict.fregre.lm <- function(object, newdata = NULL, new.scalar = NULL, ...) {
  if (is.null(newdata)) {
    return(object$fitted.values)
  }
  if (!inherits(newdata, "fdata")) {
    newdata <- fdata(newdata, argvals = object$fdataobj$argvals)
  }

  # Project new data into FPC space
  X_centered <- sweep(newdata$data, 2, object$.fpca_mean)
  scores <- X_centered %*% object$.fpca_rotation

  ncomp <- object$ncomp
  coefs <- object$coefficients
  # coefficients[1] is the intercept; FPC coefficients start at index 2
  fpc_coefs <- coefs[seq_len(ncomp) + 1L]
  y_pred <- rep(object$intercept, nrow(newdata$data))
  y_pred <- y_pred + as.vector(scores[, seq_len(ncomp), drop = FALSE] %*%
                                 fpc_coefs)

  if (!is.null(new.scalar) && length(object$gamma) > 0) {
    new.scalar <- as.matrix(new.scalar)
    y_pred <- y_pred + as.vector(new.scalar %*% object$gamma)
  }

  y_pred
}

#' Bootstrap Confidence Intervals for Functional Coefficient
#'
#' Computes pointwise and simultaneous bootstrap confidence intervals for
#' the functional coefficient beta(t) in a functional linear or logistic model.
#'
#' @param model A fitted object of class 'fregre.lm' or 'fregre.logistic'.
#' @param n.boot Number of bootstrap replicates (default 200).
#' @param alpha Significance level (default 0.05 for 95 percent CI).
#' @param seed Random seed for reproducibility.
#'
#' @return An object of class 'fregre.bootstrap.ci' with components:
#' \itemize{
#'   \item \code{lower}, \code{upper}: Pointwise CI bounds.
#'   \item \code{center}: Original beta(t) estimate.
#'   \item \code{sim.lower}, \code{sim.upper}: Simultaneous CI bands.
#'   \item \code{n.boot.success}: Number of successful bootstrap replicates.
#'   \item \code{argvals}: Grid points.
#'   \item \code{alpha}: Significance level used.
#'   \item \code{model.class}: Class of the original model.
#' }
#'
#' @examples
#' \dontrun{
#' fit <- fregre.lm(fdataobj, y, ncomp = 3)
#' ci <- fregre.bootstrap.ci(fit, n.boot = 200)
#' plot(ci)
#' }
#'
#' @export
fregre.bootstrap.ci <- function(model, n.boot = 200, alpha = 0.05, seed = 42) {
  is_logistic <- inherits(model, "fregre.logistic")
  is_lm <- inherits(model, "fregre.lm")
  if (!is_logistic && !is_lm) {
    stop("model must be of class 'fregre.lm' or 'fregre.logistic'")
  }

  fdataobj <- model$fdataobj
  y <- model$y
  ncomp <- model$ncomp
  data_mat <- fdataobj$data

  sc <- model$scalar.covariates
  if (!is.null(sc)) {
    sc <- as.matrix(sc)
  }

  if (is_logistic) {
    max_iter <- if (!is.null(model$iterations)) model$iterations * 2L else 100L
    tol <- 1e-6
    result <- bootstrap_ci_functional_logistic_rust(
      data_mat, y, sc, as.integer(ncomp), as.integer(n.boot),
      alpha, as.integer(seed), as.integer(max_iter), tol
    )
  } else {
    result <- bootstrap_ci_fregre_lm_rust(
      data_mat, y, sc, as.integer(ncomp), as.integer(n.boot),
      alpha, as.integer(seed)
    )
  }

  if (is.null(result)) {
    stop("Bootstrap CI computation failed")
  }

  structure(
    list(
      lower = result$lower,
      upper = result$upper,
      center = result$center,
      sim.lower = result$sim_lower,
      sim.upper = result$sim_upper,
      n.boot.success = result$n_boot_success,
      argvals = fdataobj$argvals,
      alpha = alpha,
      n.boot = n.boot,
      model.class = if (is_logistic) "fregre.logistic" else "fregre.lm"
    ),
    class = "fregre.bootstrap.ci"
  )
}

#' @export
print.fregre.bootstrap.ci <- function(x, ...) {
  cat("Bootstrap Confidence Intervals for beta(t)\n")
  cat("  Model type:", x$model.class, "\n")
  cat("  Alpha:", x$alpha, "\n")
  cat("  Bootstrap replicates:", x$n.boot.success, "/", x$n.boot, "succeeded\n")
  cat("  Grid points:", length(x$argvals), "\n")
  invisible(x)
}

#' @export
plot.fregre.bootstrap.ci <- function(x, simultaneous = FALSE, ...) {
  argvals <- x$argvals
  center <- x$center
  if (simultaneous) {
    lower <- x$sim.lower
    upper <- x$sim.upper
    band_label <- "Simultaneous"
  } else {
    lower <- x$lower
    upper <- x$upper
    band_label <- "Pointwise"
  }

  df <- data.frame(
    t = argvals,
    beta = center,
    lower = lower,
    upper = upper
  )

  ci_pct <- paste0(round((1 - x$alpha) * 100), "%")
  ggplot2::ggplot(df, ggplot2::aes(x = .data$t)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
                         fill = "steelblue", alpha = 0.3) +
    ggplot2::geom_line(ggplot2::aes(y = .data$beta), linewidth = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::labs(
      x = "t", y = expression(beta(t)),
      title = paste(band_label, ci_pct, "CI for", expression(beta(t)))
    ) +
    ggplot2::theme_minimal()
}

#' Random Projection Statistic
#'
#' Computes Cramer-von Mises and Kolmogorov-Smirnov statistics for
#' random projection goodness-of-fit testing.
#'
#' @param proj.x.ord Integer vector of projection orderings.
#' @param residuals Numeric vector of residuals.
#' @param n.proj Number of random projections.
#'
#' @return A list with \code{cvm} and \code{ks} vectors of test statistics.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 10
#' n_proj <- 3
#' proj <- unlist(replicate(n_proj, sample(n), simplify = FALSE))
#' rp.stat(proj, rnorm(n), n.proj = n_proj)
#' }
#'
#' @export
rp.stat <- function(proj.x.ord, residuals, n.proj) {
  result <- rp_stat(as.integer(proj.x.ord), as.double(residuals),
                    as.integer(n.proj))
  list(cvm = result$cvm, ks = result$ks)
}

# =============================================================================
# Robust Regression
# =============================================================================

#' L1 (Least Absolute Deviation) Functional Regression
#'
#' Fits a functional linear model using L1 loss (median regression),
#' which is robust to outliers in the response.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param y Response vector.
#' @param scalar.covariates Optional matrix of scalar covariates.
#' @param ncomp Number of FPC components (default 3).
#'
#' @return An object of class 'fregre.robust' with components:
#' \describe{
#'   \item{intercept}{Intercept}
#'   \item{beta.t}{Functional coefficient beta(t)}
#'   \item{fitted.values}{Fitted response values}
#'   \item{residuals}{Residuals}
#'   \item{coefficients}{Regression coefficients}
#'   \item{ncomp}{Number of FPC components used}
#'   \item{iterations}{Number of IRLS iterations}
#'   \item{converged}{Whether the algorithm converged}
#'   \item{weights}{Final IRLS weights}
#'   \item{r.squared}{R-squared statistic}
#'   \item{method}{Either "l1" or "huber"}
#' }
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' fit <- fregre.l1(fd, y, ncomp = 3)
#' }
#'
#' @export
fregre.l1 <- function(fdataobj, y, scalar.covariates = NULL, ncomp = 3) {
  if (!inherits(fdataobj, "fdata")) stop("fdataobj must be of class 'fdata'")

  sc <- if (!is.null(scalar.covariates)) as.matrix(scalar.covariates) else NULL

  result <- .Call("wrap__fregre_l1_rust",
                  fdataobj$data, as.numeric(y), sc, as.integer(ncomp))

  if (is.null(result)) stop("fregre.l1 failed")

  out <- list(
    intercept = result$intercept,
    beta.t = result$beta_t,
    fitted.values = result$fitted_values,
    residuals = result$residuals,
    coefficients = result$coefficients,
    ncomp = result$ncomp,
    iterations = result$iterations,
    converged = result$converged,
    weights = result$weights,
    r.squared = result$r_squared,
    method = "l1",
    .fpca_mean = result$fpca_mean,
    .fpca_rotation = matrix(result$fpca_rotation_data,
                            nrow = result$fpca_rotation_nrow,
                            ncol = result$fpca_rotation_ncol),
    .fpca_scores = matrix(result$fpca_scores_data,
                          nrow = result$fpca_scores_nrow,
                          ncol = result$fpca_scores_ncol),
    fdataobj = fdataobj,
    y = y
  )
  class(out) <- "fregre.robust"
  out
}

#' Huber M-Estimation Functional Regression
#'
#' Fits a functional linear model using Huber's M-estimation,
#' which provides a smooth compromise between least squares and L1.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param y Response vector.
#' @param scalar.covariates Optional matrix of scalar covariates.
#' @param ncomp Number of FPC components (default 3).
#' @param k Huber tuning parameter (default 1.345 for 95 percent efficiency).
#'
#' @return An object of class 'fregre.robust'. See \code{\link{fregre.l1}}.
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' fit <- fregre.huber(fd, y, ncomp = 3)
#' }
#'
#' @export
fregre.huber <- function(fdataobj, y, scalar.covariates = NULL, ncomp = 3, k = 1.345) {
  if (!inherits(fdataobj, "fdata")) stop("fdataobj must be of class 'fdata'")

  sc <- if (!is.null(scalar.covariates)) as.matrix(scalar.covariates) else NULL

  result <- .Call("wrap__fregre_huber_rust",
                  fdataobj$data, as.numeric(y), sc, as.integer(ncomp),
                  as.numeric(k))

  if (is.null(result)) stop("fregre.huber failed")

  out <- list(
    intercept = result$intercept,
    beta.t = result$beta_t,
    fitted.values = result$fitted_values,
    residuals = result$residuals,
    coefficients = result$coefficients,
    ncomp = result$ncomp,
    iterations = result$iterations,
    converged = result$converged,
    weights = result$weights,
    r.squared = result$r_squared,
    method = "huber",
    k = k,
    .fpca_mean = result$fpca_mean,
    .fpca_rotation = matrix(result$fpca_rotation_data,
                            nrow = result$fpca_rotation_nrow,
                            ncol = result$fpca_rotation_ncol),
    .fpca_scores = matrix(result$fpca_scores_data,
                          nrow = result$fpca_scores_nrow,
                          ncol = result$fpca_scores_ncol),
    fdataobj = fdataobj,
    y = y
  )
  class(out) <- "fregre.robust"
  out
}

#' Predict from Robust Functional Regression
#'
#' @param object A fitted 'fregre.robust' model.
#' @param newdata An object of class 'fdata' (new curves).
#' @param scalar.covariates Optional new scalar covariates.
#' @param ... Additional arguments (ignored).
#'
#' @return Predicted response values.
#'
#' @examples
#' \donttest{
#' fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
#' y <- rnorm(50)
#' fit <- fregre.l1(fd, y, ncomp = 3)
#' pred <- predict(fit, fd[1:5, ])
#' }
#'
#' @export
predict.fregre.robust <- function(object, newdata, scalar.covariates = NULL, ...) {
  if (!inherits(newdata, "fdata")) stop("newdata must be of class 'fdata'")

  sc <- if (!is.null(scalar.covariates)) as.matrix(scalar.covariates) else NULL

  .Call("wrap__predict_fregre_robust_rust",
        object$intercept, as.numeric(object$coefficients),
        as.numeric(object$.fpca_mean),
        as.numeric(object$.fpca_rotation),
        as.integer(nrow(object$.fpca_rotation)),
        as.integer(ncol(object$.fpca_rotation)),
        as.integer(object$ncomp), newdata$data, sc)
}

#' @export
print.fregre.robust <- function(x, ...) {
  cat("Robust Functional Linear Model (", x$method, ")\n", sep = "")
  cat("  Components:", x$ncomp, "\n")
  cat("  R-squared:", round(x$r.squared, 4), "\n")
  cat("  IRLS iterations:", x$iterations, "\n")
  cat("  Converged:", x$converged, "\n")
  invisible(x)
}
