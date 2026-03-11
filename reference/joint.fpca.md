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
