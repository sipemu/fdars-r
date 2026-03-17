# Elastic Clustering for Functional Data

Functions for clustering functional data using elastic (Fisher-Rao)
distances, including k-means and hierarchical clustering in the elastic
metric. Elastic K-Means Clustering

## Usage

``` r
elastic.kmeans(
  fdataobj,
  k,
  max.iter = 100,
  tol = 1e-04,
  karcher.max.iter = 20,
  karcher.tol = 1e-04,
  lambda = 0,
  seed = 42
)
```

## Arguments

- fdataobj:

  An object of class `fdata`.

- k:

  Integer: number of clusters.

- max.iter:

  Maximum number of k-means iterations (default 100).

- tol:

  Convergence tolerance on assignment changes (default 1e-4).

- karcher.max.iter:

  Maximum iterations for Karcher mean within each cluster (default 20).

- karcher.tol:

  Convergence tolerance for Karcher mean (default 1e-4).

- lambda:

  Regularization parameter for elastic alignment (default 0).

- seed:

  Random seed for initialization (default 42).

## Value

An object of class `elastic.kmeans` with components:

- labels:

  Integer vector of cluster assignments (1-indexed)

- centers:

  List of cluster centers, each a list with `mean`, `mean_srsf`,
  `gammas`, `aligned_data`, `n_iter`, `converged`

- within.distances:

  Within-cluster distance for each observation

- total.within.distance:

  Total within-cluster distance

- n.iter:

  Number of k-means iterations

- converged:

  Logical: did the algorithm converge?

- fdataobj:

  Original fdata object

- k:

  Number of clusters

- call:

  The matched call

## Details

Partition functional data into k clusters using elastic (Fisher-Rao)
distances. Cluster centers are computed as Karcher means within each
cluster. The algorithm alternates between assignment and center update
steps until convergence.

## References

Srivastava, A. and Klassen, E. (2016). *Functional and Shape Data
Analysis*. Springer.

Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
functional data using phase and amplitude separation. *Computational
Statistics & Data Analysis*, 61:50–66.

## See also

[`elastic.hclust`](https://sipemu.github.io/fdars-r/reference/elastic.hclust.md)
for hierarchical clustering,
[`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md)
for Karcher means,
[`elastic.distance`](https://sipemu.github.io/fdars-r/reference/elastic.distance.md)
for elastic distances

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 30)
X <- matrix(0, 12, 30)
for (i in 1:6)  X[i, ] <- sin(2 * pi * t) + rnorm(30, 0, 0.1)
for (i in 7:12) X[i, ] <- cos(2 * pi * t) + rnorm(30, 0, 0.1)
fd <- fdata(X, argvals = t)
cl <- elastic.kmeans(fd, k = 2)
cl
#> Elastic K-Means Clustering
#>   Curves: 12 x 30 grid points
#>   Clusters: 2 
#>   Iterations: 1 
#>   Converged: TRUE 
#>   Total within-cluster distance: 7.581 
#>   Cluster sizes: 6, 6 
# }
```
