# Joint (Amplitude + Phase) FPCA

Performs joint FPCA combining amplitude and phase components with a
balancing parameter.

## Usage

``` r
joint.fpca(karcher, ncomp = 3)
```

## Arguments

- karcher:

  A Karcher mean result from
  [`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md).

- ncomp:

  Number of principal components (default 3).

## Value

An object of class 'joint.fpca' with components:

- scores:

  Matrix of joint FPC scores (n x ncomp)

- eigenvalues:

  Eigenvalues

- cumulative.variance:

  Cumulative proportion of variance

- balance.c:

  Balancing parameter between amplitude and phase

- vert.component:

  Vertical (amplitude) eigenvector components

- horiz.component:

  Horizontal (phase) eigenvector components

- karcher:

  The Karcher mean used

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
km <- karcher.mean(fd)
j <- joint.fpca(km, ncomp = 2)
j
#> Joint (Amplitude + Phase) FPCA
#>   Components: 2 
#>   Balance parameter: 0.000331 
#>   Cumulative variance:  59 %, 100 % 
# }
```
