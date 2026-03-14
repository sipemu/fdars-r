# SPM Contribution Diagnostics

When an alarm is triggered, identifies which functional variables or
components contribute most to the elevated T-squared or SPE statistic.
Useful for root-cause analysis.

## Usage

``` r
spm.contributions(
  scores,
  eigenvalues = NULL,
  grid.sizes = NULL,
  type = c("t2", "spe"),
  standardized.vars = NULL,
  reconstructed.vars = NULL,
  argvals.list = NULL
)
```

## Arguments

- scores:

  Score matrix (n x ncomp).

- eigenvalues:

  Eigenvalues vector (length ncomp). Required for T-squared.

- grid.sizes:

  Integer vector of grid sizes per variable. For univariate SPM, use a
  single value equal to ncomp.

- type:

  Character; either `"t2"` for T-squared contributions or `"spe"` for
  SPE contributions. Default `"t2"`.

- standardized.vars:

  For SPE contributions: list of standardized data matrices.

- reconstructed.vars:

  For SPE contributions: list of reconstructed data matrices.

- argvals.list:

  For SPE contributions: list of argvals vectors.

## Value

A matrix (n x p) of per-variable contributions.

## Examples

``` r
# \donttest{
set.seed(1)
scores <- matrix(rnorm(10 * 3), 10, 3)
eigenvalues <- c(5, 2, 1)
# Per-component contributions (univariate, 1 variable with 3 grid points)
contrib <- spm.contributions(scores, eigenvalues, grid.sizes = 3L)
dim(contrib)
#> [1] 10  1
# }
```
