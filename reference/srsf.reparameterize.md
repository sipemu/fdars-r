# Apply Warping Function to a Curve

Reparameterizes a curve f by composing it with a warping function gamma:
f(gamma(t)).

## Usage

``` r
srsf.reparameterize(fdataobj, gamma)
```

## Arguments

- fdataobj:

  An fdata object (single curve, first row used).

- gamma:

  Numeric vector of warping function values on the same grid.

## Value

An fdata object containing the reparameterized curve.
