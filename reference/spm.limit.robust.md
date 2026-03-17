# Robust Control Limits via Alternative Methods

Computes T-squared or SPE control limits using methods beyond the
standard parametric approach: empirical quantile, bootstrap, or kernel
density estimation. Useful when the parametric chi-squared assumption is
questionable.

## Usage

``` r
spm.limit.robust(
  values,
  ncomp = NULL,
  alpha = 0.05,
  type = c("t2", "spe"),
  method = c("parametric", "empirical", "bootstrap", "kde"),
  n.bootstrap = 1000,
  seed = 42
)
```

## Arguments

- values:

  Numeric vector of Phase I statistic values (T-squared or SPE).

- ncomp:

  Number of components (required for T-squared; ignored for SPE).

- alpha:

  Significance level (default 0.05).

- type:

  Which statistic: `"t2"` or `"spe"`.

- method:

  Estimation method: `"parametric"`, `"empirical"`, `"bootstrap"`, or
  `"kde"`.

- n.bootstrap:

  Number of bootstrap resamples (default 1000; only used when
  `method = "bootstrap"`).

- seed:

  Random seed for bootstrap (default 42).

## Value

A list with components:

- ucl:

  Upper control limit

- alpha:

  Significance level used

- description:

  Character string describing the method

- method:

  The method used

## See also

[`spm.phase1`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md)
for standard control limits

## Examples

``` r
# \donttest{
set.seed(1)
t2_values <- rchisq(100, df = 3)
spm.limit.robust(t2_values, ncomp = 3, type = "t2", method = "empirical")
#> $ucl
#> [1] 5.587954
#> 
#> $alpha
#> [1] 0.05
#> 
#> $description
#> [1] "T2 empirical quantile, alpha=0.05"
#> 
#> $method
#> [1] "empirical"
#> 
spm.limit.robust(t2_values, ncomp = 3, type = "t2", method = "bootstrap")
#> $ucl
#> [1] 5.587954
#> 
#> $alpha
#> [1] 0.05
#> 
#> $description
#> [1] "T2 bootstrap (1000 resamples), alpha=0.05"
#> 
#> $method
#> [1] "bootstrap"
#> 
# }
```
