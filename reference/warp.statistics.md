# Warping Function Statistics

Compute summary statistics for a collection of warping functions from a
Karcher mean computation, including the mean warp, variance, standard
deviation bands, and geodesic distances.

## Usage

``` r
warp.statistics(karcher.result)
```

## Arguments

- karcher.result:

  An object of class `karcher.mean`.

## Value

An object of class `warp.statistics` with components:

- mean:

  Mean warping function

- variance:

  Pointwise variance of warps

- std.dev:

  Pointwise standard deviation of warps

- lower.band:

  Lower confidence band (mean - 2 SD)

- upper.band:

  Upper confidence band (mean + 2 SD)

- karcher.mean.warp:

  Karcher mean of warps in geodesic space

- geodesic.distances:

  Geodesic distances of each warp from identity

- call:

  The matched call

## References

Srivastava, A. and Klassen, E. (2016). *Functional and Shape Data
Analysis*. Springer.

## See also

[`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md),
[`phase.boxplot`](https://sipemu.github.io/fdars-r/reference/phase.boxplot.md)

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 50)
X <- matrix(0, 15, 50)
for (i in 1:15) X[i, ] <- sin(2 * pi * (t - i / 60))
fd <- fdata(X, argvals = t)
km <- karcher.mean(fd, max.iter = 5)
ws <- warp.statistics(km)
ws
#> Warping Function Statistics
#>   Grid points: 50 
#>   Curves: 15 
#>   Mean geodesic distance: 0.1262 
#>   Max geodesic distance: 0.2318 
# }
```
