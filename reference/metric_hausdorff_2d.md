# Compute Hausdorff distance for 2D functional data (surfaces)

For surfaces, each sample is treated as a point cloud in 3D space: (s_i,
t_j, f(s_i, t_j)) : for all grid points

## Usage

``` r
metric_hausdorff_2d(fdata, argvals_s, argvals_t)
```

## Details

The Hausdorff distance measures how far apart two such point clouds are.
