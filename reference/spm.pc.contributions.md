# Per-PC Contributions to T-Squared Statistic

Computes the contribution of each principal component to the Hotelling
T-squared statistic for each observation. Useful for diagnosing which
components drive an alarm.

## Usage

``` r
spm.pc.contributions(scores, eigenvalues)
```

## Arguments

- scores:

  Score matrix (n x ncomp), e.g., from `spm.monitor()$scores`.

- eigenvalues:

  Numeric vector of eigenvalues (length ncomp).

## Value

A matrix (n x ncomp) of per-component T-squared contributions. Row sums
equal the T-squared statistics.

## See also

[`spm.contributions`](https://sipemu.github.io/fdars-r/reference/spm.contributions.md)
for per-variable contributions

## Examples

``` r
# \donttest{
set.seed(1)
scores <- matrix(rnorm(10 * 3), 10, 3)
eigenvalues <- c(5, 2, 1)
contrib <- spm.pc.contributions(scores, eigenvalues)
dim(contrib)  # 10 x 3
#> [1] 10  3
# }
```
