# Bayesian Pairwise Curve Alignment

Align two curves using a Bayesian framework with MCMC sampling of the
warping function posterior. This provides uncertainty quantification for
the alignment via posterior samples of the warping function.

## Usage

``` r
bayesian.align.pair(
  f1,
  f2,
  argvals,
  n.samples = 2000,
  burn.in = 500,
  step.size = 0.1,
  proposal.variance = 1,
  seed = 42
)
```

## Arguments

- f1:

  Numeric vector of the first curve's values.

- f2:

  Numeric vector of the second curve's values.

- argvals:

  Numeric vector of the evaluation grid (common to both curves).

- n.samples:

  Number of MCMC samples to draw (default 2000).

- burn.in:

  Number of initial samples to discard as burn-in (default 500).

- step.size:

  Step size for the MCMC proposal (default 0.1).

- proposal.variance:

  Variance of the MCMC proposal distribution (default 1.0).

- seed:

  Random seed for reproducibility (default 42).

## Value

A list with components:

- gamma.mean:

  numeric vector of the posterior mean warping function

- gamma.samples:

  matrix of posterior warping function samples (samples x grid)

- f2.aligned:

  numeric vector of the aligned second curve

- acceptance.rate:

  the MCMC acceptance rate

## References

Cheng, W., Dryden, I.L., and Huang, X. (2016). Bayesian registration of
functions and curves. *Bayesian Analysis*, 11(2):447–475.

## See also

[`elastic.align`](https://sipemu.github.io/fdars-r/reference/elastic.align.md)

## Examples

``` r
# \donttest{
t <- seq(0, 1, length.out = 100)
f1 <- sin(2 * pi * t)
f2 <- sin(2 * pi * (t^1.3))
res <- bayesian.align.pair(f1, f2, t, n.samples = 500, burn.in = 100)
# }
```
