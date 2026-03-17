# Elastic Hierarchical Clustering

Perform agglomerative hierarchical clustering on functional data using
pairwise elastic (Fisher-Rao) distances. Supports single, complete, and
average linkage.

## Usage

``` r
elastic.hclust(
  fdataobj,
  method = c("complete", "single", "average"),
  lambda = 0
)
```

## Arguments

- fdataobj:

  An object of class `fdata`.

- method:

  Character: linkage method. One of `"complete"` (default), `"single"`,
  or `"average"`.

- lambda:

  Regularization parameter for elastic alignment (default 0).

## Value

An object of class `elastic.hclust` with components:

- merges:

  A 3-column matrix of merge steps (i, j, distance), with 1-indexed
  cluster labels

- distance.matrix:

  The pairwise elastic distance matrix

- method:

  The linkage method used

- fdataobj:

  Original fdata object

- call:

  The matched call

## References

Srivastava, A. and Klassen, E. (2016). *Functional and Shape Data
Analysis*. Springer.

## See also

[`elastic.cutree`](https://sipemu.github.io/fdars-r/reference/elastic.cutree.md)
for cutting the dendrogram,
[`elastic.kmeans`](https://sipemu.github.io/fdars-r/reference/elastic.kmeans.md)
for k-means clustering,
[`elastic.distance`](https://sipemu.github.io/fdars-r/reference/elastic.distance.md)
for the distance metric

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 30)
X <- matrix(0, 10, 30)
for (i in 1:5)  X[i, ] <- sin(2 * pi * t) + rnorm(30, 0, 0.1)
for (i in 6:10) X[i, ] <- cos(2 * pi * t) + rnorm(30, 0, 0.1)
fd <- fdata(X, argvals = t)
hc <- elastic.hclust(fd, method = "complete")
hc
#> Elastic Hierarchical Clustering
#>   Curves: 10 x 30 grid points
#>   Method: complete 
#>   Merge steps: 9 
# }
```
