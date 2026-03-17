# Cut Elastic Dendrogram

Cut an elastic hierarchical clustering dendrogram to produce a specified
number of clusters.

## Usage

``` r
elastic.cutree(tree, k)
```

## Arguments

- tree:

  An object of class `elastic.hclust`.

- k:

  Integer: desired number of clusters.

## Value

An integer vector of cluster labels (1-indexed), one per curve.

## See also

[`elastic.hclust`](https://sipemu.github.io/fdars-r/reference/elastic.hclust.md)

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
labels <- elastic.cutree(hc, k = 2)
table(labels)
#> labels
#> 1 2 
#> 5 5 
# }
```
