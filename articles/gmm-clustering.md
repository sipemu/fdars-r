# Model-Based Clustering with Gaussian Mixtures

## Introduction

Gaussian Mixture Models (GMMs) provide model-based clustering with
several advantages over k-means:

- **Soft assignment**: Each curve gets a probability of belonging to
  each cluster, not just a hard label.
- **Automatic model selection**: BIC or ICL criteria select the number
  of clusters.
- **Flexible covariance**: Full, diagonal, or spherical covariance
  structures adapt to cluster shapes.

**fdars** implements GMM clustering through
[`cluster.gmm()`](https://sipemu.github.io/fdars-r/reference/cluster.gmm.md),
which projects functional data onto a B-spline basis before fitting the
mixture model.

``` r
library(fdars)
#> 
#> Attaching package: 'fdars'
#> The following objects are masked from 'package:stats':
#> 
#>     cov, decompose, deriv, median, sd, var
#> The following object is masked from 'package:base':
#> 
#>     norm
library(ggplot2)
theme_set(theme_minimal())
```

## Simulated Data

``` r
set.seed(42)
n_per_cluster <- 25
m <- 100
t_grid <- seq(0, 1, length.out = m)

# Cluster 1: Sine curves
X1 <- matrix(0, n_per_cluster, m)
for (i in 1:n_per_cluster) {
  X1[i, ] <- sin(2 * pi * t_grid) + rnorm(m, sd = 0.15)
}

# Cluster 2: Cosine curves
X2 <- matrix(0, n_per_cluster, m)
for (i in 1:n_per_cluster) {
  X2[i, ] <- cos(2 * pi * t_grid) + rnorm(m, sd = 0.15)
}

# Cluster 3: Linear curves
X3 <- matrix(0, n_per_cluster, m)
for (i in 1:n_per_cluster) {
  X3[i, ] <- 2 * t_grid - 1 + rnorm(m, sd = 0.15)
}

X <- rbind(X1, X2, X3)
true_clusters <- rep(1:3, each = n_per_cluster)

fd <- fdata(X, argvals = t_grid)
plot(fd)
```

![](gmm-clustering_files/figure-html/data-1.png)

## Fitting a GMM

[`cluster.gmm()`](https://sipemu.github.io/fdars-r/reference/cluster.gmm.md)
searches over a range of k values and selects the best model using BIC
(default) or ICL.

``` r
gmm_fit <- cluster.gmm(fd, k.range = 2:6, seed = 123)
print(gmm_fit)
#> Gaussian Mixture Model Clustering
#> ==================================
#>   Number of clusters: 6 
#>   Number of observations: 75 
#>   Criterion: BIC 
#>   BIC: 322.95 
#>   ICL: 322.95 
#>   Converged: TRUE 
#>   Covariance type: full 
#>   Basis: bspline ( 10 functions )
#> 
#> Cluster sizes:
#> 
#>  1  2  3  4  5  6 
#> 16 14 21  9  9  6
```

## Model Selection: BIC Curve

The BIC values across candidate k help visualize model selection:

``` r
plot(gmm_fit)
```

![](gmm-clustering_files/figure-html/bic-plot-1.png)

## Soft vs Hard Assignment

Unlike k-means, GMM provides membership probabilities:

``` r
# Show membership for first 6 curves
head(round(gmm_fit$membership, 3))
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    0    0    0    0    0
#> [2,]    0    0    1    0    0    0
#> [3,]    0    0    0    1    0    0
#> [4,]    0    0    0    0    0    1
#> [5,]    0    0    1    0    0    0
#> [6,]    0    0    1    0    0    0

# Maximum membership gives hard assignment
cat("Hard clusters:", head(gmm_fit$cluster), "\n")
#> Hard clusters: 1 3 4 6 3 3
```

## Comparison with K-Means

``` r
kmeans_fit <- cluster.kmeans(fd, ncl = 3, seed = 123)

cat("K-means WCSS:", round(kmeans_fit$tot.withinss, 2), "\n")
#> K-means WCSS: 1.63
cat("GMM BIC:", round(gmm_fit$bic, 2), "\n")
#> GMM BIC: 322.95
cat("GMM converged:", gmm_fit$converged, "\n")
#> GMM converged: TRUE
```

## Covariance Types

[`cluster.gmm()`](https://sipemu.github.io/fdars-r/reference/cluster.gmm.md)
supports different covariance structures:

- `"full"` (default): Each cluster has its own full covariance matrix.
- `"diagonal"`: Diagonal covariance (features are independent within
  clusters).

``` r
gmm_diag <- cluster.gmm(fd, k.range = 3, cov.type = "diagonal", seed = 123)
cat("Full GMM BIC:", round(gmm_fit$bic, 2), "\n")
#> Full GMM BIC: 322.95
cat("Diagonal GMM BIC:", round(gmm_diag$bic, 2), "\n")
#> Diagonal GMM BIC: 1729.75
```

## Predicting New Observations

Assign new curves to clusters using the fitted model:

``` r
# Generate new curves
set.seed(99)
X_new <- matrix(0, 3, m)
X_new[1, ] <- sin(2 * pi * t_grid) + rnorm(m, sd = 0.1)   # Should be cluster 1
X_new[2, ] <- cos(2 * pi * t_grid) + rnorm(m, sd = 0.1)   # Should be cluster 2
X_new[3, ] <- 2 * t_grid - 1 + rnorm(m, sd = 0.1)          # Should be cluster 3

fd_new <- fdata(X_new, argvals = t_grid)
pred <- predict(gmm_fit, fd_new)

cat("Predicted clusters:", pred$cluster, "\n")
#> Predicted clusters: 2 2 2
cat("Membership probabilities:\n")
#> Membership probabilities:
print(round(pred$membership, 3))
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    0    1    0    0    0    0
#> [2,]    0    1    0    0    0    0
#> [3,]    0    1    0    0    0    0
```

## BIC vs ICL

ICL (Integrated Completed Likelihood) penalizes cluster overlap more
strongly than BIC:

``` r
gmm_icl <- cluster.gmm(fd, k.range = 2:6, criterion = "icl", seed = 123)
cat("BIC selected k:", gmm_fit$k, "\n")
#> BIC selected k: 6
cat("ICL selected k:", gmm_icl$k, "\n")
#> ICL selected k: 6
```

## See Also

- [`vignette("articles/clustering")`](https://sipemu.github.io/fdars-r/articles/clustering.md)
  for k-means and fuzzy c-means
- [`vignette("articles/functional-classification")`](https://sipemu.github.io/fdars-r/articles/functional-classification.md)
  for supervised classification with
  [`fclassif()`](https://sipemu.github.io/fdars-r/reference/fclassif.md)
