# Pairwise Alignment Consistency

Measure the consistency of pairwise elastic alignments by checking the
triangle closure property. For each triplet (i, j, k), consistency
checks whether aligning i to k directly gives the same result as
aligning i to j then j to k.

## Usage

``` r
alignment.pairwise.consistency(fdataobj, lambda = 0, max.triplets = 0)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- lambda:

  Penalty weight on warp deviation from identity. Default 0.

- max.triplets:

  Maximum number of random triplets to evaluate. Default 0 (use all
  triplets).

## Value

A scalar in \[0, 1\] measuring alignment consistency.

## Details

A value near 1 indicates high consistency (transitive alignments); a
value near 0 suggests the data contains distinct subgroups or the
alignment is sensitive to the target choice.

## Examples

``` r
# \donttest{
t <- seq(0, 1, length.out = 100)
X <- matrix(0, 10, 100)
for (i in 1:10) X[i, ] <- sin(2*pi*(t - i/50))
fd <- fdata(X, argvals = t)
alignment.pairwise.consistency(fd, lambda = 0)
#> [1] 0.004180664
# }
```
