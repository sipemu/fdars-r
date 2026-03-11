# Prototype and Criticism

Finds representative (prototype) and atypical (criticism) observations
using the Maximum Mean Discrepancy (MMD) framework on FPCA scores.

## Usage

``` r
fregre.prototype(model, ncomp, n.prototypes = 5, n.criticisms = 5)
```

## Arguments

- model:

  A fitted `fregre.lm` or `fregre.logistic` model (uses its FPCA
  scores).

- ncomp:

  Number of components to use.

- n.prototypes:

  Number of prototypes (default 5).

- n.criticisms:

  Number of criticisms (default 5).

## Value

A list with prototypes (indices), prototype_witness, criticisms
(indices), criticism_witness, and bandwidth.
