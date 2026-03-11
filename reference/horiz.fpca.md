# Horizontal (Phase) FPCA

Performs FPCA on the phase (warping) component after elastic alignment.
Captures timing variability independent of shape.

## Usage

``` r
horiz.fpca(karcher, ncomp = 3)
```

## Arguments

- karcher:

  A Karcher mean result from
  [`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md).

- ncomp:

  Number of principal components (default 3).

## Value

An object of class 'horiz.fpca' with components:

- scores:

  Matrix of FPC scores (n x ncomp)

- eigenfunctions.psi:

  Eigenfunctions in tangent space

- eigenfunctions.gam:

  Eigenfunctions as warp functions

- eigenvalues:

  Eigenvalues

- cumulative.variance:

  Cumulative proportion of variance

- mean.psi:

  Mean tangent vector

- shooting.vectors:

  Shooting vectors matrix

- karcher:

  The Karcher mean used
