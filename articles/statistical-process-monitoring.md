# Statistical Process Monitoring for Functional Data

## Introduction

Statistical Process Monitoring (SPM) extends classical control charts to
functional data. When the “quality characteristic” is an entire curve
(e.g., a spectral profile, a temperature trajectory, or a process
signature), univariate control charts on summary statistics lose
information. Functional SPM monitors the full curve structure.

![SPM Pipeline: Phase I estimation and Phase II
monitoring](../reference/figures/spm-overview.svg)

The workflow has two phases:

1.  **Phase I** — Estimate the in-control process using historical
    “good” data: perform FPCA, compute Hotelling $T^{2}$ and Squared
    Prediction Error (SPE) statistics, and set control limits.

2.  **Phase II** — Monitor new observations by projecting onto the Phase
    I principal components and comparing statistics to the control
    limits. Alarms fire when a statistic exceeds its Upper Control Limit
    (UCL).

### When to Use Functional SPM

| Scenario                           | Method                                                                                                                                                    |
|------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------|
| Single functional quality variable | [`spm.phase1()`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md) + [`spm.monitor()`](https://sipemu.github.io/fdars-r/reference/spm.monitor.md) |
| Multiple functional variables      | `mf.spm.phase1()` + `mf.spm.monitor()`                                                                                                                    |
| Detecting gradual drifts           | [`spm.ewma()`](https://sipemu.github.io/fdars-r/reference/spm.ewma.md)                                                                                    |
| Process with known covariates      | `frcc()` + [`frcc.monitor()`](https://sipemu.github.io/fdars-r/reference/frcc.monitor.md)                                                                 |
| Diagnosing alarm root cause        | [`spm.contributions()`](https://sipemu.github.io/fdars-r/reference/spm.contributions.md)                                                                  |

## Simulating Process Data

We simulate semiconductor wafer thickness profiles: each curve is a
spectral measurement across the wafer surface. Phase I uses 200
in-control profiles; Phase II introduces a shift at observation 30.

``` r
library(fdars)
set.seed(42)

m <- 50  # grid points
t_grid <- seq(0, 1, length.out = m)

# In-control process: smooth curves with random variation
n_ic <- 200  # Phase I observations
X_ic <- matrix(0, n_ic, m)
for (i in seq_len(n_ic)) {
  scores <- rnorm(3, sd = c(2, 1, 0.5))
  X_ic[i, ] <- 5 + scores[1] * sin(2 * pi * t_grid) +
                    scores[2] * cos(4 * pi * t_grid) +
                    scores[3] * sin(6 * pi * t_grid) +
                    rnorm(m, sd = 0.3)
}

# Phase II: 50 observations, shift at observation 30
n_new <- 50
X_new <- matrix(0, n_new, m)
for (i in seq_len(n_new)) {
  scores <- rnorm(3, sd = c(2, 1, 0.5))
  shift <- if (i >= 30) 1.5 * sin(2 * pi * t_grid) else 0
  X_new[i, ] <- 5 + scores[1] * sin(2 * pi * t_grid) +
                     scores[2] * cos(4 * pi * t_grid) +
                     scores[3] * sin(6 * pi * t_grid) +
                     shift + rnorm(m, sd = 0.3)
}

fd_ic <- fdata(X_ic, argvals = t_grid)
fd_new <- fdata(X_new, argvals = t_grid)
```

## Phase I: Estimating Control Limits

[`spm.phase1()`](https://sipemu.github.io/fdars-r/reference/spm.phase1.md)
performs FPCA on the in-control data, computes $T^{2}$ and SPE
statistics for each training observation, and estimates control limits
at significance level $\alpha$.

``` r
chart <- spm.phase1(fd_ic, ncomp = 3, alpha = 0.05, seed = 42)
print(chart)
```

The returned `spm.chart` object contains:

- **T2 control limit** — based on the $F$-distribution with `ncomp` and
  `n - ncomp` degrees of freedom
- **SPE control limit** — estimated from the empirical distribution of
  Phase I SPE values
- **Phase I statistics** — $T^{2}$ and SPE for each training observation
  (useful for checking that Phase I data is truly in-control)

## Phase II: Online Monitoring

Project new observations onto the Phase I principal components and
compare to control limits:

``` r
monitor <- spm.monitor(chart, fd_new)
print(monitor)

# Which observations triggered alarms?
alarm_idx <- which(monitor$t2.alarm | monitor$spe.alarm)
cat("Alarms at observations:", alarm_idx, "\n")
```

The monitoring result contains:

| Field       | Description                                       |
|-------------|---------------------------------------------------|
| `t2`        | Hotelling $T^{2}$ values for each new observation |
| `spe`       | SPE values for each new observation               |
| `t2.alarm`  | Logical: did $T^{2}$ exceed the UCL?              |
| `spe.alarm` | Logical: did SPE exceed the UCL?                  |
| `scores`    | PC scores for the new observations                |

## EWMA Monitoring

The Exponentially Weighted Moving Average (EWMA) chart smooths the
monitoring statistics over time, making it more sensitive to small,
persistent shifts than the standard Shewhart-type chart:

``` r
ewma_result <- spm.ewma(chart, fd_new, lambda = 0.2, alpha = 0.05)

cat("EWMA T2 alarms:", sum(ewma_result$t2.alarm), "\n")
cat("EWMA SPE alarms:", sum(ewma_result$spe.alarm), "\n")
```

The smoothing parameter $\lambda \in (0,1\rbrack$ controls the memory: -
$\lambda = 1$ gives the standard Shewhart chart (no smoothing) -
$\lambda = 0.2$ (default) gives moderate smoothing - Smaller $\lambda$
detects smaller shifts but responds more slowly

## Multivariate FPCA

When multiple functional variables are measured simultaneously (e.g.,
temperature, pressure, and flow rate profiles), Multivariate FPCA
extracts joint modes of variation across all variables:

``` r
# Simulate two correlated functional variables
fd_var1 <- fdata(X_ic, argvals = t_grid)
fd_var2 <- fdata(X_ic + matrix(rnorm(n_ic * m, sd = 0.5), n_ic, m),
                 argvals = t_grid)

mfpca_result <- mfpca(list(fd_var1, fd_var2), ncomp = 3)
cat("Eigenvalues:", round(mfpca_result$eigenvalues, 3), "\n")
```

The result contains joint scores and per-variable eigenfunctions, which
can be used for multivariate SPM via `mf.spm.phase1()`.

## Contribution Diagnostics

When an alarm fires, contribution analysis identifies *which* region of
the curve or *which* variable caused the alarm:

``` r
# T2 contributions show which PC scores are responsible
t2_contrib <- spm.contributions(monitor$scores,
                                chart$eigenvalues,
                                chart$grid.sizes,
                                type = "t2")

# Each row = observation, each column = variable/region contribution
cat("Contribution of PC1 to observation 35:", round(t2_contrib[35, 1], 3), "\n")
```

## Comparison of Methods

| Method   | Detects                                                  | Best for                             |
|----------|----------------------------------------------------------|--------------------------------------|
| **T2**   | Mean shifts in the PC score space                        | Location changes in dominant modes   |
| **SPE**  | Changes not captured by the PCs                          | New types of variation, sensor drift |
| **EWMA** | Small persistent shifts                                  | Gradual process degradation          |
| **FRCC** | Shifts in the relationship between curves and covariates | Regression-based monitoring          |

## See Also

- [`vignette("articles/fpca")`](https://sipemu.github.io/fdars-r/articles/fpca.md)
  — Functional PCA fundamentals
- [`vignette("articles/streaming-depth")`](https://sipemu.github.io/fdars-r/articles/streaming-depth.md)
  — Online depth-based monitoring
- [`vignette("articles/outlier-detection")`](https://sipemu.github.io/fdars-r/articles/outlier-detection.md)
  — Depth-based outlier detection

## References

- Colosimo, B.M. and Pacella, M. (2010). A comparison study of control
  charts for statistical monitoring of functional data. *International
  Journal of Production Research*, 48(6), 1575–1601.
- Paynabar, K., Jin, J., and Pacella, M. (2013). Monitoring and
  diagnosis of multichannel nonlinear profiles using uncorrelated
  multilinear principal component analysis. *IIE Transactions*, 45(11),
  1235–1247.
- Grasso, M., Colosimo, B.M., and Pacella, M. (2014). Profile monitoring
  via sensor fusion: the use of PCA methods for multi-channel data.
  *International Journal of Production Research*, 52(20), 6110–6135.
