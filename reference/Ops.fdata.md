# Arithmetic Operations for Functional Data

Perform elementwise arithmetic on `fdata` objects. Supports addition,
subtraction, multiplication, division, and exponentiation between two
`fdata` objects or between an `fdata` object and a numeric scalar.

## Usage

``` r
# S3 method for class 'fdata'
Ops(e1, e2)
```

## Arguments

- e1, e2:

  Objects of class `fdata` or numeric scalars. At least one must be of
  class `fdata`.

## Value

An `fdata` object with the result.

## Examples

``` r
fd1 <- fdata(matrix(1:20, 4, 5))
fd2 <- fdata(matrix(21:40, 4, 5))
fd1 + fd2
#> Functional data object
#>   Type: 1D (curve) 
#>   Number of observations: 4 
#>   Number of points: 5 
#>   Range: 1 - 5 
fd1 * 2
#> Functional data object
#>   Type: 1D (curve) 
#>   Number of observations: 4 
#>   Number of points: 5 
#>   Range: 1 - 5 
3 - fd1
#> Functional data object
#>   Type: 1D (curve) 
#>   Number of observations: 4 
#>   Number of points: 5 
#>   Range: 1 - 5 
```
