# DFBETAS and DFFITS

Computes DFBETAS (influence on each coefficient) and DFFITS (influence
on fitted values) for each observation.

## Usage

``` r
fregre.dfbetas(model, data)
```

## Arguments

- model:

  A fitted `fregre.lm` model.

- data:

  An fdata object (the training data).

## Value

A list with dfbetas (matrix), dffits, studentized_residuals,
dfbetas_cutoff, and dffits_cutoff.

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), nrow = 50), argvals = seq(0, 1, length.out = 10))
y <- rnorm(50)
fit <- fregre.lm(fd, y, ncomp = 3)
result <- fregre.dfbetas(fit, fd)
# }
```
