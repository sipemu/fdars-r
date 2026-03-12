# Beta Decomposition

Decomposes the beta coefficient function into per-FPC contributions.

## Usage

``` r
fregre.beta.decomp(model)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model.

## Value

A list with components (per-FPC beta contributions), coefficients, and
variance_proportion.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
result <- fregre.beta.decomp(fit)
# }
```
