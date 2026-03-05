# Compute Fraiman-Muniz depth

Uses the FM1 formula from fda.usc: d = 1 - \|0.5 - Fn(x)\| With
scale=TRUE (default): d = (1 - \|0.5 - Fn(x)\| - 0.5) \* 2 = 2 \*
min(Fn(x), 1-Fn(x))

## Usage

``` r
depth_fm_1d(fdataobj, fdataori, `_trim`, scale)
```
