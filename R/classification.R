#' Supervised Classification of Functional Data
#'
#' Functions for classifying functional observations into discrete groups.

#' Functional Classification
#'
#' Classifies functional data using one of several methods: LDA, QDA, kNN,
#' kernel, or DD-plot (depth-vs-depth).
#'
#' @param fdataobj An object of class 'fdata'.
#' @param y Integer vector of class labels (1-indexed).
#' @param method Classification method: "lda", "qda", "knn", "kernel", or "dd".
#' @param covariates Optional matrix of scalar covariates.
#' @param ncomp Number of FPC components (default 3). Used by lda, qda, knn.
#' @param ... Additional arguments:
#'   \describe{
#'     \item{k}{Number of neighbors for kNN (default 5).}
#'     \item{h.func}{Bandwidth for functional kernel (default auto).}
#'     \item{h.scalar}{Bandwidth for scalar kernel (default auto).}
#'   }
#'
#' @return An object of class 'fclassif' with components:
#'   \item{predicted}{Predicted class labels}
#'   \item{probabilities}{Posterior probabilities (if available)}
#'   \item{accuracy}{Training accuracy}
#'   \item{confusion}{Confusion matrix}
#'   \item{method}{Method used}
#'   \item{ncomp}{Number of FPC components}
#'
#' @export
fclassif <- function(fdataobj, y, method = c("lda", "qda", "knn", "kernel", "dd"),
                     covariates = NULL, ncomp = 3, ...) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  method <- match.arg(method)
  n <- nrow(fdataobj$data)
  y <- as.integer(y)

  if (length(y) != n) {
    stop("Length of y must equal number of curves")
  }

  cov_mat <- if (!is.null(covariates)) as.matrix(covariates) else NULL
  dots <- list(...)

  result <- switch(method,
    "lda" = .Call("wrap__fclassif_lda_rust", fdataobj$data, y, cov_mat,
                  as.integer(ncomp)),
    "qda" = .Call("wrap__fclassif_qda_rust", fdataobj$data, y, cov_mat,
                  as.integer(ncomp)),
    "knn" = {
      k <- dots$k %||% 5L
      .Call("wrap__fclassif_knn_rust", fdataobj$data, y, cov_mat,
            as.integer(ncomp), as.integer(k))
    },
    "kernel" = {
      h.func <- dots$h.func %||% 0.0
      h.scalar <- dots$h.scalar %||% 0.0
      .Call("wrap__fclassif_kernel_rust", fdataobj$data,
            as.numeric(fdataobj$argvals), y, cov_mat,
            as.numeric(h.func), as.numeric(h.scalar))
    },
    "dd" = .Call("wrap__fclassif_dd_rust", fdataobj$data, y, cov_mat)
  )

  if (is.null(result)) {
    stop("fclassif (", method, ") failed: check data dimensions and labels")
  }

  structure(
    list(
      predicted = result$predicted,
      probabilities = if (!is.null(result$probabilities) &&
                          !inherits(result$probabilities, "NULL")) {
        result$probabilities
      } else {
        NULL
      },
      accuracy = result$accuracy,
      confusion = result$confusion,
      method = method,
      ncomp = result$ncomp,
      n.classes = result$n_classes,
      fdataobj = fdataobj,
      y = y,
      covariates = covariates,
      call = match.call()
    ),
    class = "fclassif"
  )
}

#' Cross-Validated Functional Classification
#'
#' Evaluates classification error rate using k-fold cross-validation.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param y Integer vector of class labels.
#' @param method Classification method (default "lda").
#' @param covariates Optional scalar covariates matrix.
#' @param ncomp Number of FPC components (default 3).
#' @param nfold Number of CV folds (default 10).
#' @param seed Random seed for fold assignment.
#'
#' @return An object of class 'fclassif.cv' with components:
#'   \item{error.rate}{Mean error rate across folds}
#'   \item{fold.errors}{Per-fold error rates}
#'   \item{best.ncomp}{Best ncomp if tuned}
#'
#' @export
fclassif.cv <- function(fdataobj, y, method = "lda", covariates = NULL,
                         ncomp = 3, nfold = 10, seed = NULL) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  y <- as.integer(y)
  cov_mat <- if (!is.null(covariates)) as.matrix(covariates) else NULL
  seed_val <- if (!is.null(seed)) as.integer(seed) else 42L

  result <- .Call("wrap__fclassif_cv_rust", fdataobj$data,
                  as.numeric(fdataobj$argvals), y, cov_mat,
                  as.character(method), as.integer(ncomp),
                  as.integer(nfold), as.integer(seed_val))

  if (is.null(result)) {
    stop("fclassif.cv failed")
  }

  structure(
    list(
      error.rate = result$error_rate,
      fold.errors = result$fold_errors,
      best.ncomp = result$best_ncomp,
      method = method,
      nfold = nfold
    ),
    class = "fclassif.cv"
  )
}

#' Print method for fclassif objects
#' @param x A fclassif object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object \code{x}.
#' @export
print.fclassif <- function(x, ...) {
  cat("Functional Classification\n")
  cat("=========================\n")
  cat("  Method:", toupper(x$method), "\n")
  cat("  Number of observations:", length(x$y), "\n")
  cat("  Number of classes:", x$n.classes, "\n")
  if (!is.null(x$ncomp) && x$ncomp > 0) {
    cat("  FPC components:", x$ncomp, "\n")
  }
  cat("  Training accuracy:", round(x$accuracy, 4), "\n")
  cat("\nConfusion matrix:\n")
  print(x$confusion)
  invisible(x)
}

#' Print method for fclassif.cv objects
#' @param x A fclassif.cv object.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object \code{x}.
#' @export
print.fclassif.cv <- function(x, ...) {
  cat("Cross-Validated Functional Classification\n")
  cat("==========================================\n")
  cat("  Method:", toupper(x$method), "\n")
  cat("  Folds:", x$nfold, "\n")
  cat("  Error rate:", round(x$error.rate, 4), "\n")
  cat("  Best ncomp:", x$best.ncomp, "\n")
  invisible(x)
}

#' Plot method for fclassif objects
#'
#' Plots the confusion matrix as a heatmap.
#'
#' @param x A fclassif object.
#' @param ... Additional arguments (ignored).
#' @return A ggplot object.
#' @export
plot.fclassif <- function(x, ...) {
  conf <- x$confusion
  nc <- nrow(conf)

  df <- data.frame(
    true = rep(seq_len(nc), each = nc),
    predicted = rep(seq_len(nc), nc),
    count = as.vector(t(conf))
  )

  ggplot2::ggplot(df, ggplot2::aes(x = factor(.data$predicted),
                                    y = factor(.data$true),
                                    fill = .data$count)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = .data$count), color = "white", size = 5) +
    ggplot2::scale_fill_gradient(low = "steelblue", high = "darkred") +
    ggplot2::labs(x = "Predicted", y = "True",
                  title = paste0("Confusion Matrix (", toupper(x$method),
                                 ", acc=", round(x$accuracy, 3), ")"),
                  fill = "Count") +
    ggplot2::theme(legend.position = "none")
}
