# Modified Epigraph Index (MEI) for 1D functional data MEI measures the proportion of time a curve is below other curves MEI(x_i) = (1/n) \* sum_j (1/m) \* sum_t I(x_i(t) \< x_j(t)) + 0.5\*I(x_i(t) = x_j(t))

Modified Epigraph Index (MEI) for 1D functional data MEI measures the
proportion of time a curve is below other curves MEI(x_i) = (1/n) \*
sum_j (1/m) \* sum_t I(x_i(t) \< x_j(t)) + 0.5\*I(x_i(t) = x_j(t))

## Usage

``` r
depth_mei_1d(fdataobj, fdataori)
```
