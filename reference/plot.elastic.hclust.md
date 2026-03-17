# Plot Elastic Hierarchical Clustering Result

Displays a dendrogram of the hierarchical clustering. The merge
information is converted to an `hclust` object for plotting.

## Usage

``` r
# S3 method for class 'elastic.hclust'
plot(x, ...)
```

## Arguments

- x:

  An object of class `elastic.hclust`.

- ...:

  Additional arguments passed to `plot.hclust`.

## Value

Invisible NULL. Called for its side effect of producing a plot.
