# Unified Outlier Summary

Produces a ranked data frame from any outlier detection result, showing
the outlierness score, p-value, and outlier status for each curve.

## Usage

``` r
outlier_summary(x)
```

## Arguments

- x:

  An outlier detection result: an object of class `"outliers.fdata"`,
  `"magnitudeshape"`, or `"outliergram"`.

## Value

A `data.frame` with columns:

- index:

  Curve index (1-based)

- outlierness:

  Continuous outlierness score (higher = more outlying)

- p.value:

  P-value for the outlier test

- is_outlier:

  Logical indicating whether the curve was flagged

Rows are sorted by outlierness score in descending order.

## Examples

``` r
set.seed(42)
t <- seq(0, 1, length.out = 50)
X <- matrix(0, 30, 50)
for (i in 1:28) X[i, ] <- sin(2*pi*t) + rnorm(50, sd = 0.2)
X[29, ] <- sin(2*pi*t) + 3
X[30, ] <- -sin(2*pi*t)
fd <- fdata(X, argvals = t)

out <- outliers.depth.pond(fd, nb = 100)
outlier_summary(out)
#>    index outlierness    p.value is_outlier
#> 29    29  0.96774194 0.03225806       TRUE
#> 30    30  0.93548387 0.06451613       TRUE
#> 20    20  0.90322581 0.09677419      FALSE
#> 11    11  0.87096774 0.12903226      FALSE
#> 1      1  0.83870968 0.16129032      FALSE
#> 21    21  0.80645161 0.19354839      FALSE
#> 17    17  0.77419355 0.22580645      FALSE
#> 25    25  0.74193548 0.25806452      FALSE
#> 6      6  0.70967742 0.29032258      FALSE
#> 13    13  0.67741935 0.32258065      FALSE
#> 19    19  0.64516129 0.35483871      FALSE
#> 24    24  0.61290323 0.38709677      FALSE
#> 18    18  0.58064516 0.41935484      FALSE
#> 10    10  0.54838710 0.45161290      FALSE
#> 9      9  0.51612903 0.48387097      FALSE
#> 12    12  0.48387097 0.51612903      FALSE
#> 27    27  0.45161290 0.54838710      FALSE
#> 22    22  0.41935484 0.58064516      FALSE
#> 5      5  0.38709677 0.61290323      FALSE
#> 26    26  0.35483871 0.64516129      FALSE
#> 8      8  0.32258065 0.67741935      FALSE
#> 2      2  0.29032258 0.70967742      FALSE
#> 3      3  0.25806452 0.74193548      FALSE
#> 23    23  0.22580645 0.77419355      FALSE
#> 16    16  0.19354839 0.80645161      FALSE
#> 14    14  0.16129032 0.83870968      FALSE
#> 4      4  0.12903226 0.87096774      FALSE
#> 15    15  0.09677419 0.90322581      FALSE
#> 28    28  0.06451613 0.93548387      FALSE
#> 7      7  0.03225806 0.96774194      FALSE
```
