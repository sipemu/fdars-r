# Elastic Distance Matrix

Compute the pairwise elastic (Fisher-Rao) distance matrix for functional
data.

## Usage

``` r
elastic.distance(fdataobj, fdataref = NULL, ...)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- fdataref:

  Optional reference 'fdata'. If NULL, computes self-distances.

- ...:

  Additional arguments (ignored).

## Value

A distance matrix (numeric matrix).

## References

Srivastava, A., Klassen, E., Joshi, S.H., and Jermyn, I.H. (2011). Shape
analysis of elastic curves in Euclidean spaces. *IEEE Transactions on
Pattern Analysis and Machine Intelligence*, 33(7):1415–1428.

Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
functional data using phase and amplitude separation. *Computational
Statistics & Data Analysis*, 61:50–66.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(100), 10, 10), argvals = seq(0, 1, length.out = 10))
D <- elastic.distance(fd)
# }
```
