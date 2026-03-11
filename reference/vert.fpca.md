# Elastic FPCA

Amplitude, phase, and joint FPCA after elastic alignment using the
Square Root Velocity Framework (SRVF). Vertical (Amplitude) FPCA

## Usage

``` r
vert.fpca(karcher, ncomp = 3)
```

## Arguments

- karcher:

  A Karcher mean result from
  [`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md).

- ncomp:

  Number of principal components (default 3).

## Value

An object of class 'vert.fpca' with components:

- scores:

  Matrix of FPC scores (n x ncomp)

- eigenfunctions.q:

  SRVF eigenfunctions

- eigenfunctions.f:

  Original-space eigenfunctions

- eigenvalues:

  Eigenvalues (variance explained by each component)

- cumulative.variance:

  Cumulative proportion of variance

- mean.q:

  Mean SRVF

- karcher:

  The Karcher mean used

## Details

Performs FPCA on the amplitude component of curves after elastic
alignment. Captures shape variability independent of timing/phase.
