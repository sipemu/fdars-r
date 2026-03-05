# Nadaraya-Watson smoother matrix S_ij = K((t_i - t_j)/h) \* w_j / sum_k(K((t_i - t_k)/h) \* w_k)

Nadaraya-Watson smoother matrix S_ij = K((t_i - t_j)/h) \* w_j /
sum_k(K((t_i - t_k)/h) \* w_k)

## Usage

``` r
s_nw(argvals, h, kernel_type, weights, cv)
```
