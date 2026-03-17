# Phase Boxplot for Warping Functions

Compute a functional boxplot for warping functions, identifying the
median, central 50% region, whiskers, and outliers based on functional
depth.

## Usage

``` r
phase.boxplot(karcher.result, factor = 1.5)
```

## Arguments

- karcher.result:

  An object of class `karcher.mean`.

- factor:

  Whisker extension factor, analogous to the IQR multiplier in a
  classical boxplot (default 1.5).

## Value

An object of class `phase.boxplot` with components:

- median:

  Median warping function

- median.index:

  Index of the median warp in the input

- central.lower:

  Lower boundary of the central 50% region

- central.upper:

  Upper boundary of the central 50% region

- whisker.lower:

  Lower whisker

- whisker.upper:

  Upper whisker

- depths:

  Functional depth values for each warp

- outlier.indices:

  Indices of outlying warps

- factor:

  The whisker extension factor used

- argvals:

  Grid points

- call:

  The matched call

## References

Sun, Y. and Genton, M.G. (2011). Functional boxplots. *Journal of
Computational and Graphical Statistics*, 20(2):316–334.

## See also

[`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md),
[`warp.statistics`](https://sipemu.github.io/fdars-r/reference/warp.statistics.md)

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 50)
X <- matrix(0, 15, 50)
for (i in 1:15) X[i, ] <- sin(2 * pi * (t - i / 60))
fd <- fdata(X, argvals = t)
km <- karcher.mean(fd, max.iter = 5)
pb <- phase.boxplot(km)
pb
#> Phase Boxplot (Warping Functions)
#>   Grid points: 50 
#>   Curves: 15 
#>   Median index: 8 
#>   Outliers: 2 
#>   Factor: 1.5 
# }
```
