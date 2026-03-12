# Bootstrap Confidence Intervals for Functional Coefficient

Computes pointwise and simultaneous bootstrap confidence intervals for
the functional coefficient beta(t) in a functional linear or logistic
model.

## Usage

``` r
fregre.bootstrap.ci(model, n.boot = 200, alpha = 0.05, seed = 42)
```

## Arguments

- model:

  A fitted object of class 'fregre.lm' or 'fregre.logistic'.

- n.boot:

  Number of bootstrap replicates (default 200).

- alpha:

  Significance level (default 0.05 for 95 percent CI).

- seed:

  Random seed for reproducibility.

## Value

An object of class 'fregre.bootstrap.ci' with components:

- `lower`, `upper`: Pointwise CI bounds.

- `center`: Original beta(t) estimate.

- `sim.lower`, `sim.upper`: Simultaneous CI bands.

- `n.boot.success`: Number of successful bootstrap replicates.

- `argvals`: Grid points.

- `alpha`: Significance level used.

- `model.class`: Class of the original model.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fregre.lm(fdataobj, y, ncomp = 3)
ci <- fregre.bootstrap.ci(fit, n.boot = 200)
plot(ci)
} # }
```
