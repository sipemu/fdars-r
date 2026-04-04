# Elastic Depth for Functional Data

Compute elastic depth scores that decompose into amplitude and phase
components. Elastic depth measures the centrality of each curve relative
to the sample, accounting for both vertical variation (amplitude) and
horizontal variation (phase/warping).

## Usage

``` r
elastic.depth(fdataobj, lambda = 0)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- lambda:

  Regularisation parameter controlling warping smoothness (default 0).

## Value

An object of class 'elastic.depth' with components:

- amplitude_depth:

  numeric vector of amplitude depth scores

- phase_depth:

  numeric vector of phase depth scores

- combined_depth:

  numeric vector of combined depth scores

## References

Lopez-Pintado, S. and Romo, J. (2009). On the concept of depth for
functional data. *Journal of the American Statistical Association*,
104(486):718–734.

Srivastava, A. and Klassen, E. (2016). *Functional and Shape Data
Analysis*. Springer.

## See also

[`elastic.outlier.detection`](https://sipemu.github.io/fdars-r/reference/elastic.outlier.detection.md)

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 50)
X <- matrix(0, 15, 50)
for (i in 1:15) X[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
fd <- fdata(X, argvals = t)
ed <- elastic.depth(fd)
# }
```
