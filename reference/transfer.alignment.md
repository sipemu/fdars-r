# Transfer Alignment Between Functional Samples

Transfer the alignment learned from a source sample of curves to a
target sample. The Karcher mean and warping structure estimated from the
source are used to guide the alignment of the target curves, which is
useful when the target sample is small or noisy.

## Usage

``` r
transfer.alignment(source, target, lambda = 0, max.iter = 15, tol = 1e-04)
```

## Arguments

- source:

  An object of class 'fdata' (the reference sample).

- target:

  An object of class 'fdata' (the sample to align).

- lambda:

  Regularisation parameter controlling warping smoothness (default 0).

- max.iter:

  Maximum number of iterations (default 15).

- tol:

  Convergence tolerance (default 1e-4).

## Value

A list with components:

- aligned:

  fdata of the aligned target curves

- gammas:

  fdata of warping functions applied to the target

- source.mean:

  fdata of the source Karcher mean

- distances:

  numeric vector of elastic distances

## See also

[`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md),
[`elastic.align`](https://sipemu.github.io/fdars-r/reference/elastic.align.md)

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 50)
X1 <- matrix(0, 20, 50)
for (i in 1:20) X1[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
X2 <- matrix(0, 5, 50)
for (i in 1:5) X2[i, ] <- sin(2 * pi * t + runif(1, -0.5, 0.5))
source <- fdata(X1, argvals = t)
target <- fdata(X2, argvals = t)
res <- transfer.alignment(source, target)
# }
```
