# Tolerance Bands for Functional Data

Compute tolerance bands that are expected to contain a given fraction of
individual curves in the population. Functional Tolerance Band

## Usage

``` r
tolerance.band(
  fdataobj,
  method = c("fpca", "conformal", "scb", "exponential", "elastic", "phase",
    "elastic.config"),
  coverage = 0.95,
  ncomp = 3,
  nb = 500,
  band.type = c("pointwise", "simultaneous"),
  cal.fraction = 0.2,
  score.type = c("supnorm", "l2"),
  bandwidth = NULL,
  confidence = NULL,
  multiplier = c("gaussian", "rademacher"),
  family = c("gaussian", "binomial", "poisson"),
  max.iter = 10,
  ncomp.phase = 3,
  tol = 1e-04,
  seed = NULL
)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- method:

  Method to use. One of "fpca" (default), "conformal", "scb",
  "exponential", "elastic", "phase", or "elastic.config".

- coverage:

  Target coverage probability (default 0.95).

- ncomp:

  Number of FPCA components (default 3). Used by "fpca", "exponential",
  "elastic", "phase", and "elastic.config" (as amplitude components for
  "elastic.config").

- nb:

  Number of bootstrap replicates (default 500). Used by "fpca", "scb",
  "exponential", "elastic", "phase", and "elastic.config".

- band.type:

  "pointwise" (default) or "simultaneous". Used by "fpca", "elastic",
  "phase", and "elastic.config".

- cal.fraction:

  Calibration fraction for conformal method (default 0.2).

- score.type:

  Nonconformity score: "supnorm" (default) or "l2". Used by "conformal".

- bandwidth:

  Kernel bandwidth for SCB Degras method. If NULL, a default is
  computed.

- confidence:

  Confidence level for SCB method (default is `coverage`).

- multiplier:

  Multiplier distribution: "gaussian" (default) or "rademacher". Used by
  "scb".

- family:

  Exponential family: "gaussian" (default), "binomial", or "poisson".
  Used by "exponential".

- max.iter:

  Maximum iterations for elastic method Karcher mean (default 10).

- ncomp.phase:

  Number of FPCA components for the phase band (default 3). Used by
  "elastic.config".

- tol:

  Convergence tolerance for Karcher mean (default 1e-4). Used by
  "elastic.config".

- seed:

  Random seed for reproducibility (default NULL).

## Value

An object of class 'tolerance.band' with components:

- lower:

  numeric vector of lower bounds

- upper:

  numeric vector of upper bounds

- center:

  numeric vector of center function

- half_width:

  numeric vector of half-widths

- method:

  the method used

- coverage:

  the target coverage

- argvals:

  evaluation points

- fdataobj:

  the original fdata input

The "phase" method additionally returns `gamma.lower`, `gamma.upper`,
and `gamma.center` (warping function bounds), with
`lower`/`upper`/`center`/`half_width` from the tangent-space band.

The "elastic.config" method additionally returns `phase.lower`,
`phase.upper`, and `phase.center` (warping function bounds), with the
primary `lower`/`upper`/`center`/`half_width` from the amplitude band.

Returns NULL with a warning if computation fails.

## Details

Compute a tolerance band for functional data using one of several
methods.

Available methods:

- fpca:

  FPCA + bootstrap tolerance band. Reconstructs curves from PC scores
  and uses bootstrap to estimate pointwise or simultaneous quantiles.

- conformal:

  Distribution-free conformal prediction band. Splits data into training
  and calibration sets.

- scb:

  Simultaneous confidence band for the mean (Degras method). Uses
  multiplier bootstrap for critical values.

- exponential:

  Tolerance band for exponential family functional data. Applies link
  function transformation.

- elastic:

  Tolerance band in elastic (aligned) space. First computes Karcher
  mean, then applies FPCA band on aligned data.

- phase:

  Phase tolerance band for warping function variation. Computes Karcher
  mean, extracts warping functions, and builds a tolerance band in the
  tangent (shooting vector) space. Returns both the warping-function
  bounds and the tangent-space band.

- elastic.config:

  Joint amplitude and phase tolerance band with full configuration
  control. Separately controls the number of FPCA components for
  amplitude (`ncomp`) and phase (`ncomp.phase`), plus convergence
  tolerance (`tol`).

## References

Rathnayake, L.N. and Cuevas, A. (2016). Tolerance bands for functional
data. *Technometrics*, 58(3):326–334.

Lei, J. and Wasserman, L. (2014). Distribution-free prediction bands for
non-parametric regression. *Journal of the Royal Statistical Society:
Series B*, 76(1):71–96.

Degras, D. (2011). Simultaneous confidence bands for nonparametric
regression with functional data. *Statistica Sinica*, 21(4):1735–1765.

## Examples

``` r
fd <- fdata(matrix(rnorm(500), 50, 10), argvals = seq(0, 1, length.out = 10))
# \donttest{
band <- tolerance.band(fd, method = "fpca", nb = 100)
# }
```
