#' Elastic Clustering for Functional Data
#'
#' Functions for clustering functional data using elastic (Fisher-Rao) distances,
#' including k-means and hierarchical clustering in the elastic metric.

# =============================================================================
# Elastic K-Means
# =============================================================================

#' Elastic K-Means Clustering
#'
#' Partition functional data into k clusters using elastic (Fisher-Rao)
#' distances. Cluster centers are computed as Karcher means within each
#' cluster. The algorithm alternates between assignment and center update
#' steps until convergence.
#'
#' @param fdataobj An object of class \code{fdata}.
#' @param k Integer: number of clusters.
#' @param max.iter Maximum number of k-means iterations (default 100).
#' @param tol Convergence tolerance on assignment changes (default 1e-4).
#' @param karcher.max.iter Maximum iterations for Karcher mean within each
#'   cluster (default 20).
#' @param karcher.tol Convergence tolerance for Karcher mean (default 1e-4).
#' @param lambda Regularization parameter for elastic alignment (default 0).
#' @param seed Random seed for initialization (default 42).
#'
#' @return An object of class \code{elastic.kmeans} with components:
#' \describe{
#'   \item{labels}{Integer vector of cluster assignments (1-indexed)}
#'   \item{centers}{List of cluster centers, each a list with
#'     \code{mean}, \code{mean_srsf}, \code{gammas}, \code{aligned_data},
#'     \code{n_iter}, \code{converged}}
#'   \item{within.distances}{Within-cluster distance for each observation}
#'   \item{total.within.distance}{Total within-cluster distance}
#'   \item{n.iter}{Number of k-means iterations}
#'   \item{converged}{Logical: did the algorithm converge?}
#'   \item{fdataobj}{Original fdata object}
#'   \item{k}{Number of clusters}
#'   \item{call}{The matched call}
#' }
#'
#' @references
#' Srivastava, A. and Klassen, E. (2016). \emph{Functional and Shape Data
#' Analysis}. Springer.
#'
#' Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
#' functional data using phase and amplitude separation.
#' \emph{Computational Statistics & Data Analysis}, 61:50--66.
#'
#' @seealso \code{\link{elastic.hclust}} for hierarchical clustering,
#'   \code{\link{karcher.mean}} for Karcher means,
#'   \code{\link{elastic.distance}} for elastic distances
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 30)
#' X <- matrix(0, 12, 30)
#' for (i in 1:6)  X[i, ] <- sin(2 * pi * t) + rnorm(30, 0, 0.1)
#' for (i in 7:12) X[i, ] <- cos(2 * pi * t) + rnorm(30, 0, 0.1)
#' fd <- fdata(X, argvals = t)
#' cl <- elastic.kmeans(fd, k = 2)
#' cl
#' }
elastic.kmeans <- function(fdataobj, k, max.iter = 100, tol = 1e-4,
                            karcher.max.iter = 20, karcher.tol = 1e-4,
                            lambda = 0, seed = 42) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }

  k <- as.integer(k)
  n <- nrow(fdataobj$data)

  if (k < 1 || k > n) {
    stop("k must be between 1 and the number of curves")
  }

  argvals <- as.numeric(fdataobj$argvals)

  result <- .Call("wrap__elastic_kmeans_rust",
                  fdataobj$data, argvals,
                  k,
                  as.integer(max.iter), as.numeric(tol),
                  as.integer(karcher.max.iter), as.numeric(karcher.tol),
                  as.numeric(lambda), as.integer(seed))

  if (is.null(result)) {
    stop("elastic.kmeans failed: check data dimensions and parameters")
  }

  structure(
    list(
      labels = result$labels,
      centers = result$centers,
      within.distances = result$within_distances,
      total.within.distance = result$total_within_distance,
      n.iter = result$n_iter,
      converged = result$converged,
      fdataobj = fdataobj,
      k = k,
      call = match.call()
    ),
    class = "elastic.kmeans"
  )
}

#' Print Elastic K-Means Result
#'
#' @param x An object of class \code{elastic.kmeans}.
#' @param ... Additional arguments (ignored).
#' @export
print.elastic.kmeans <- function(x, ...) {
  n <- nrow(x$fdataobj$data)
  m <- ncol(x$fdataobj$data)
  cat("Elastic K-Means Clustering\n")
  cat("  Curves:", n, "x", m, "grid points\n")
  cat("  Clusters:", x$k, "\n")
  cat("  Iterations:", x$n.iter, "\n")
  cat("  Converged:", x$converged, "\n")
  cat("  Total within-cluster distance:",
      format(x$total.within.distance, digits = 4), "\n")
  sizes <- tabulate(x$labels, nbins = x$k)
  cat("  Cluster sizes:", paste(sizes, collapse = ", "), "\n")
  invisible(x)
}

#' Plot Elastic K-Means Result
#'
#' Displays cluster mean curves with member curves in the background,
#' one facet per cluster.
#'
#' @param x An object of class \code{elastic.kmeans}.
#' @param ... Additional arguments (ignored).
#'
#' @return A ggplot object.
#'
#' @export
plot.elastic.kmeans <- function(x, ...) {
  argvals <- as.numeric(x$fdataobj$argvals)
  m <- length(argvals)
  n <- nrow(x$fdataobj$data)
  k <- x$k

  # Build data frame of all curves with cluster labels
  df_curves <- data.frame(
    curve_id = rep(seq_len(n), each = m),
    argval = rep(argvals, n),
    value = as.vector(t(x$fdataobj$data)),
    cluster = factor(rep(x$labels, each = m))
  )

  # Build data frame of cluster means
  mean_rows <- lapply(seq_len(k), function(j) {
    center_mean <- as.numeric(x$centers[[j]]$mean)
    data.frame(
      argval = argvals,
      value = center_mean,
      cluster = factor(rep(j, m))
    )
  })
  df_means <- do.call(rbind, mean_rows)

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = df_curves,
      ggplot2::aes(x = .data$argval, y = .data$value,
                   group = .data$curve_id),
      alpha = 0.3, color = "grey60"
    ) +
    ggplot2::geom_line(
      data = df_means,
      ggplot2::aes(x = .data$argval, y = .data$value),
      color = "red", linewidth = 1.2
    ) +
    ggplot2::facet_wrap(~ .data$cluster, labeller = ggplot2::labeller(
      cluster = function(x) paste("Cluster", x)
    )) +
    ggplot2::labs(
      title = paste0("Elastic K-Means (k = ", k, ")"),
      x = "t", y = "value"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))

  p
}

# =============================================================================
# Elastic Hierarchical Clustering
# =============================================================================

#' Elastic Hierarchical Clustering
#'
#' Perform agglomerative hierarchical clustering on functional data using
#' pairwise elastic (Fisher-Rao) distances. Supports single, complete,
#' and average linkage.
#'
#' @param fdataobj An object of class \code{fdata}.
#' @param method Character: linkage method. One of \code{"complete"}
#'   (default), \code{"single"}, or \code{"average"}.
#' @param lambda Regularization parameter for elastic alignment (default 0).
#'
#' @return An object of class \code{elastic.hclust} with components:
#' \describe{
#'   \item{merges}{A 3-column matrix of merge steps (i, j, distance),
#'     with 1-indexed cluster labels}
#'   \item{distance.matrix}{The pairwise elastic distance matrix}
#'   \item{method}{The linkage method used}
#'   \item{fdataobj}{Original fdata object}
#'   \item{call}{The matched call}
#' }
#'
#' @references
#' Srivastava, A. and Klassen, E. (2016). \emph{Functional and Shape Data
#' Analysis}. Springer.
#'
#' @seealso \code{\link{elastic.cutree}} for cutting the dendrogram,
#'   \code{\link{elastic.kmeans}} for k-means clustering,
#'   \code{\link{elastic.distance}} for the distance metric
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 30)
#' X <- matrix(0, 10, 30)
#' for (i in 1:5)  X[i, ] <- sin(2 * pi * t) + rnorm(30, 0, 0.1)
#' for (i in 6:10) X[i, ] <- cos(2 * pi * t) + rnorm(30, 0, 0.1)
#' fd <- fdata(X, argvals = t)
#' hc <- elastic.hclust(fd, method = "complete")
#' hc
#' }
elastic.hclust <- function(fdataobj,
                            method = c("complete", "single", "average"),
                            lambda = 0) {
  if (!inherits(fdataobj, "fdata")) {
    stop("fdataobj must be of class 'fdata'")
  }
  method <- match.arg(method)

  n <- nrow(fdataobj$data)
  if (n < 2) {
    stop("At least 2 curves are required for hierarchical clustering")
  }

  argvals <- as.numeric(fdataobj$argvals)

  result <- .Call("wrap__elastic_hierarchical_rust",
                  fdataobj$data, argvals, method, as.numeric(lambda))

  if (is.null(result)) {
    stop("elastic.hclust failed: check data dimensions and parameters")
  }

  structure(
    list(
      merges = result$merges,
      distance.matrix = result$distance_matrix,
      method = method,
      fdataobj = fdataobj,
      call = match.call()
    ),
    class = "elastic.hclust"
  )
}

#' Print Elastic Hierarchical Clustering Result
#'
#' @param x An object of class \code{elastic.hclust}.
#' @param ... Additional arguments (ignored).
#' @export
print.elastic.hclust <- function(x, ...) {
  n <- nrow(x$fdataobj$data)
  m <- ncol(x$fdataobj$data)
  cat("Elastic Hierarchical Clustering\n")
  cat("  Curves:", n, "x", m, "grid points\n")
  cat("  Method:", x$method, "\n")
  cat("  Merge steps:", nrow(x$merges), "\n")
  invisible(x)
}

#' Plot Elastic Hierarchical Clustering Result
#'
#' Displays a dendrogram of the hierarchical clustering. The merge
#' information is converted to an \code{hclust} object for plotting.
#'
#' @param x An object of class \code{elastic.hclust}.
#' @param ... Additional arguments passed to \code{plot.hclust}.
#'
#' @return Invisible NULL. Called for its side effect of producing a plot.
#'
#' @export
plot.elastic.hclust <- function(x, ...) {
  merges <- x$merges
  n <- nrow(x$fdataobj$data)
  n_merges <- nrow(merges)

  if (n_merges == 0 || n < 2) {
    warning("No merges to plot")
    return(invisible(NULL))
  }

  # Convert merge matrix to hclust format
  # hclust$merge: negative = leaf, positive = cluster formed at that step
  # Our merges are 1-indexed (i, j, distance)
  hc_merge <- matrix(0L, n_merges, 2)
  hc_height <- numeric(n_merges)

  # Track which original indices have been merged into which step
  # Map from original 1-indexed curve to negative leaf or positive step
  cluster_map <- rep(0L, n + n_merges)
  for (i in seq_len(n)) {
    cluster_map[i] <- -i  # leaf
  }

  for (s in seq_len(n_merges)) {
    ci <- as.integer(merges[s, 1])
    cj <- as.integer(merges[s, 2])
    dist_val <- merges[s, 3]

    hc_merge[s, 1] <- cluster_map[ci]
    hc_merge[s, 2] <- cluster_map[cj]
    hc_height[s] <- dist_val

    # The new cluster is referenced by step index
    cluster_map[n + s] <- s

    # Update both merged items to point to this step
    cluster_map[ci] <- s
    cluster_map[cj] <- s
  }

  # Build hclust object
  hc <- structure(
    list(
      merge = hc_merge,
      height = hc_height,
      order = seq_len(n),
      labels = if (!is.null(rownames(x$fdataobj$data))) {
        rownames(x$fdataobj$data)
      } else {
        as.character(seq_len(n))
      },
      method = x$method,
      dist.method = "elastic"
    ),
    class = "hclust"
  )

  # Reorder for nicer dendrogram
  tryCatch({
    d <- stats::as.dist(x$distance.matrix)
    hc$order <- stats::order.dendrogram(stats::as.dendrogram(hc))
  }, error = function(e) {
    # Fall back to default order
  })

  plot(hc, main = paste("Elastic Dendrogram (", x$method, "linkage)"),
       xlab = "Curve", ylab = "Elastic Distance", ...)

  invisible(NULL)
}

# =============================================================================
# Cut Dendrogram
# =============================================================================

#' Cut Elastic Dendrogram
#'
#' Cut an elastic hierarchical clustering dendrogram to produce a specified
#' number of clusters.
#'
#' @param tree An object of class \code{elastic.hclust}.
#' @param k Integer: desired number of clusters.
#'
#' @return An integer vector of cluster labels (1-indexed), one per curve.
#'
#' @seealso \code{\link{elastic.hclust}}
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' t <- seq(0, 1, length.out = 30)
#' X <- matrix(0, 10, 30)
#' for (i in 1:5)  X[i, ] <- sin(2 * pi * t) + rnorm(30, 0, 0.1)
#' for (i in 6:10) X[i, ] <- cos(2 * pi * t) + rnorm(30, 0, 0.1)
#' fd <- fdata(X, argvals = t)
#' hc <- elastic.hclust(fd, method = "complete")
#' labels <- elastic.cutree(hc, k = 2)
#' table(labels)
#' }
elastic.cutree <- function(tree, k) {
  if (!inherits(tree, "elastic.hclust")) {
    stop("tree must be of class 'elastic.hclust'")
  }

  k <- as.integer(k)
  n <- nrow(tree$fdataobj$data)

  if (k < 1 || k > n) {
    stop("k must be between 1 and the number of curves")
  }

  result <- .Call("wrap__elastic_cut_dendrogram_rust",
                  tree$merges, tree$distance.matrix, k)

  if (is.null(result)) {
    stop("elastic.cutree failed: check tree and k")
  }

  as.integer(result)
}
