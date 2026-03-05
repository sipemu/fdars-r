#' Landmark Registration for Functional Data
#'
#' Feature-based alignment using detected landmarks (peaks, valleys,
#' zero-crossings, inflection points).

#' Detect Landmarks in Functional Data
#'
#' Detect feature landmarks (peaks, valleys, zero-crossings, inflection points)
#' in each curve of a functional data object.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param kind Character vector of landmark types to detect. One of
#'   "peak", "valley", "zero", "inflection". Default is "peak".
#' @param min.prominence Minimum prominence for peak/valley detection.
#'   Default is 0.
#'
#' @return A list of data frames (one per curve), each with columns:
#'   position, kind, value, prominence.
#'
#' @export
#' @examples
#' t <- seq(0, 1, length.out = 100)
#' X <- matrix(0, 5, 100)
#' for (i in 1:5) X[i, ] <- sin(2*pi*t + i/10)
#' fd <- fdata(X, argvals = t)
#' lms <- detect.landmarks(fd, kind = "peak")
detect.landmarks <- function(fdataobj,
                             kind = c("peak", "valley", "zero", "inflection"),
                             min.prominence = 0) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  kind <- match.arg(kind)
  n <- nrow(fdataobj$data)
  argvals <- fdataobj$argvals

  result <- vector("list", n)
  for (i in seq_len(n)) {
    res <- landmark_detect(as.numeric(fdataobj$data[i, ]),
                           as.numeric(argvals), kind,
                           as.numeric(min.prominence))
    result[[i]] <- data.frame(
      position = res$position,
      kind = res$kind,
      value = res$value,
      prominence = res$prominence,
      stringsAsFactors = FALSE
    )
  }

  result
}

#' Landmark Registration
#'
#' Register (align) functional data by detecting landmarks and warping curves
#' so that the landmarks are mapped to common target positions.
#'
#' @param fdataobj An object of class 'fdata'.
#' @param kind Character: landmark type to detect. One of
#'   "peak", "valley", "zero", "inflection". Default is "peak".
#' @param min.prominence Minimum prominence for peak/valley detection.
#'   Default is 0.
#' @param expected.count Expected number of landmarks per curve. If 0 (default),
#'   uses the minimum number detected across all curves.
#'
#' @return An object of class 'landmark.register' with components:
#' \describe{
#'   \item{registered}{fdata of registered curves}
#'   \item{gammas}{fdata of warping functions}
#'   \item{landmarks}{list of detected landmarks per curve}
#'   \item{target_landmarks}{common target landmark positions}
#'   \item{fdataobj}{the original input}
#' }
#'
#' @export
#' @examples
#' \donttest{
#' t <- seq(0, 1, length.out = 100)
#' X <- matrix(0, 10, 100)
#' for (i in 1:10) X[i, ] <- sin(2*pi*(t - i/50))
#' fd <- fdata(X, argvals = t)
#' lr <- landmark.register(fd, kind = "peak", min.prominence = 0.5)
#' }
landmark.register <- function(fdataobj,
                              kind = c("peak", "valley", "zero", "inflection"),
                              min.prominence = 0,
                              expected.count = 0) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  kind <- match.arg(kind)
  argvals <- fdataobj$argvals

  res <- landmark_register_curves(fdataobj$data, as.numeric(argvals),
                                  kind, as.numeric(min.prominence),
                                  as.integer(expected.count))

  # Detect landmarks for output
  landmarks <- detect.landmarks(fdataobj, kind = kind,
                                min.prominence = min.prominence)

  result <- list(
    registered = fdata(res$registered, argvals = argvals),
    gammas = fdata(res$gammas, argvals = argvals),
    landmarks = landmarks,
    target_landmarks = res$target_landmarks,
    fdataobj = fdataobj
  )
  class(result) <- "landmark.register"
  result
}

#' @export
print.landmark.register <- function(x, ...) {
  n <- nrow(x$fdataobj$data)
  m <- ncol(x$fdataobj$data)
  cat("Landmark Registration\n")
  cat("  Curves:", n, "x", m, "grid points\n")
  cat("  Target landmarks:", length(x$target_landmarks), "\n")
  if (length(x$target_landmarks) > 0) {
    cat("  Target positions:", paste(round(x$target_landmarks, 3), collapse = ", "), "\n")
  }
  invisible(x)
}

#' Plot Landmark Registration Results
#'
#' @param x An object of class 'landmark.register'.
#' @param type Character: "registered" (default), "original", "warps", or "both".
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object (invisibly).
#'
#' @export
plot.landmark.register <- function(x, type = c("registered", "original", "warps", "both"), ...) {
  type <- match.arg(type)

  if (type == "both") {
    p1 <- .plot_fdata_curves(x$fdataobj, title = "Original")
    p2 <- .plot_fdata_curves(x$registered, title = "Landmark Registered")
    if (requireNamespace("patchwork", quietly = TRUE)) {
      p <- p1 + p2
    } else {
      p <- p2
    }
  } else if (type == "warps") {
    p <- .plot_fdata_curves(x$gammas, title = "Warping Functions")
  } else if (type == "original") {
    p <- .plot_fdata_curves(x$fdataobj, title = "Original Curves")
  } else {
    p <- .plot_fdata_curves(x$registered, title = "Landmark Registered")
  }

  p
}
