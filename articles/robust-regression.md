# Robust Functional Regression

## Introduction

Standard least-squares regression estimates the conditional mean of the
response by minimizing
$\sum\left( y_{i} - {\widehat{y}}_{i} \right)^{2}$. This works well when
the data are clean, but even a single outlier can drastically distort
the fitted model. In functional data settings the problem is amplified:
contaminated spectra, sensor glitches, or mislabelled samples can
silently wreck an OLS fit.

**fdars** provides two robust alternatives to
[`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md): -
**[`fregre.l1()`](https://sipemu.github.io/fdars-r/reference/fregre.l1.md)**
— L1 / Least Absolute Deviation (median) regression -
**[`fregre.huber()`](https://sipemu.github.io/fdars-r/reference/fregre.huber.md)**
— Huber M-estimation, a smooth blend of L2 and L1

Both use the same FPC projection strategy as
[`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md)
but replace the OLS loss with a robust loss, solved via iteratively
reweighted least squares (IRLS). The result is an `fregre.robust` S3
object with the same interface:
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`print()`](https://rdrr.io/r/base/print.html), and access to fitted
values, residuals, and $R^{2}$.

![Robust Functional Regression: OLS vs L1 vs
Huber](../reference/figures/robust-regression-overview.svg)

### Which Method Should I Use?

| Method       | Function                                                                       | Loss                                               | Key parameter      | Best when                                   |
|--------------|--------------------------------------------------------------------------------|----------------------------------------------------|--------------------|---------------------------------------------|
| **OLS**      | [`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md)       | $\sum\left( y_{i} - {\widehat{y}}_{i} \right)^{2}$ | $K$ (# components) | Clean data; full diagnostics ecosystem      |
| **L1 (LAD)** | [`fregre.l1()`](https://sipemu.github.io/fdars-r/reference/fregre.l1.md)       | $\sum|y_{i} - {\widehat{y}}_{i}|$                  | $K$ (# components) | Heavy contamination; outliers in $y$        |
| **Huber**    | [`fregre.huber()`](https://sipemu.github.io/fdars-r/reference/fregre.huber.md) | Huber loss                                         | $K$, $k$ (tuning)  | Moderate contamination; efficiency near OLS |

**Quick decision rule:**

1.  Start with
    [`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md)
    — it is the most efficient estimator when the data are clean and has
    full downstream support (explainability, diagnostics, conformal
    prediction).
2.  If you suspect outliers in the response $y$, fit
    [`fregre.huber()`](https://sipemu.github.io/fdars-r/reference/fregre.huber.md)
    as a robust check. With the default $k = 1.345$ it retains 95%
    asymptotic efficiency under Gaussian errors.
3.  If contamination is severe (\> 10–15% of observations), use
    [`fregre.l1()`](https://sipemu.github.io/fdars-r/reference/fregre.l1.md)
    for maximum resistance.

## Simulating Clean and Contaminated Data

We simulate a near-infrared (NIR) spectroscopy scenario: 200 absorbance
spectra measured at 100 wavelengths, with a linear relationship between
the spectra and a scalar property (e.g., moisture content). We then
contaminate 10% of the responses with large measurement errors to mimic
real-world outliers.

``` r
library(fdars)
set.seed(42)

# --- Dimensions ---
n <- 200          # number of spectra
m <- 100          # number of wavelength channels
t_grid <- seq(900, 1700, length.out = m)   # NIR wavelength range (nm)

# --- Generate smooth spectra ---
# Each spectrum is a sum of three smooth components plus noise,
# mimicking overlapping absorption bands.
X <- matrix(0, n, m)
for (i in 1:n) {
  a1 <- rnorm(1, 1.0, 0.3)
  a2 <- rnorm(1, 0.5, 0.2)
  a3 <- rnorm(1, 0.3, 0.1)
  X[i, ] <- a1 * dnorm(t_grid, 1200, 80) * 200 +
            a2 * dnorm(t_grid, 1400, 60) * 200 +
            a3 * dnorm(t_grid, 1600, 50) * 200 +
            rnorm(m, sd = 0.01)
}

fd <- fdata(X, argvals = t_grid)

# --- True coefficient function ---
# Response depends on the 1200 nm and 1400 nm absorption bands.
beta_true <- 0.6 * dnorm(t_grid, 1200, 80) * 200 -
             0.4 * dnorm(t_grid, 1400, 60) * 200

# --- Clean response ---
y_clean <- numeric(n)
for (i in 1:n) {
  y_clean[i] <- sum(beta_true * X[i, ]) / m + rnorm(1, sd = 0.3)
}

# --- Contaminated response (10% outliers) ---
y_contam <- y_clean
n_outlier <- round(0.10 * n)   # 20 outliers
outlier_idx <- sample(n, n_outlier)
y_contam[outlier_idx] <- y_contam[outlier_idx] + rnorm(n_outlier, mean = 0, sd = 8)
```

The contaminated observations have errors roughly 25 times larger than
the noise standard deviation, simulating gross measurement mistakes.

## OLS Regression — Baseline

We first fit ordinary least squares via
[`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md)
on both the clean and contaminated data.

``` r
# --- Clean data ---
fit_ols_clean <- fregre.lm(fd, y_clean, ncomp = 5)
print(fit_ols_clean)

# --- Contaminated data ---
fit_ols_contam <- fregre.lm(fd, y_contam, ncomp = 5)
print(fit_ols_contam)
```

On clean data, OLS performs well. On contaminated data, the 20 outliers
pull the coefficient estimates away from their true values, inflating
the residuals for *all* observations — including the clean ones. The
$R^{2}$ drops and the RMSE increases substantially.

## L1 (LAD) Regression

[`fregre.l1()`](https://sipemu.github.io/fdars-r/reference/fregre.l1.md)
minimises the sum of absolute residuals rather than squared residuals:

$$\widehat{\mathbf{γ}} = \arg\min\limits_{\mathbf{γ}}\sum\limits_{i = 1}^{n}|y_{i} - \alpha - {\mathbf{ξ}}_{i}^{\top}{\mathbf{γ}}|$$

where ${\mathbf{ξ}}_{i}$ is the vector of FPC scores for observation
$i$. Because the absolute value penalises large residuals linearly (not
quadratically), outliers exert much less influence on the fit.

The optimisation is carried out via IRLS (iteratively reweighted least
squares), which re-weights each observation inversely proportional to
its current residual magnitude at each iteration.

``` r
# --- Clean data ---
fit_l1_clean <- fregre.l1(fd, y_clean, ncomp = 5)
print(fit_l1_clean)

# --- Contaminated data ---
fit_l1_contam <- fregre.l1(fd, y_contam, ncomp = 5)
print(fit_l1_contam)
```

On the contaminated data, the L1 fit remains close to the clean-data
fit. The IRLS weights for the 20 outlier observations will be much
smaller than those for the clean observations, effectively
down-weighting their influence.

### Examining the IRLS Weights

``` r
# Observations with small weights are effectively identified as outliers
w <- fit_l1_contam$weights
cat("Mean weight (clean obs):", round(mean(w[-outlier_idx]), 3), "\n")
cat("Mean weight (outlier obs):", round(mean(w[outlier_idx]), 3), "\n")
cat("Iterations:", fit_l1_contam$iterations, "\n")
cat("Converged:", fit_l1_contam$converged, "\n")
```

## Huber M-Estimation

Huber’s M-estimator provides a smooth compromise between L2 and L1. The
Huber loss function is:

$$\rho_{k}(r) = \begin{cases}
{\frac{1}{2}r^{2}} & {{\text{if}\mspace{6mu}}|r| \leq k} \\
{k|r| - \frac{1}{2}k^{2}} & {{\text{if}\mspace{6mu}}|r| > k}
\end{cases}$$

where $k > 0$ is a tuning parameter. For small residuals ($|r| \leq k$),
the loss is quadratic (like OLS); for large residuals ($|r| > k$), the
loss is linear (like L1). This means:

- Small residuals are treated as in standard OLS, preserving efficiency.
- Large residuals are down-weighted, providing outlier resistance.
- The default $k = 1.345$ achieves 95% asymptotic efficiency under
  Gaussian errors while still providing substantial robustness.

``` r
# --- Default k = 1.345 ---
fit_huber_contam <- fregre.huber(fd, y_contam, ncomp = 5)
print(fit_huber_contam)

# --- More aggressive robustness with smaller k ---
fit_huber_k08 <- fregre.huber(fd, y_contam, ncomp = 5, k = 0.8)
print(fit_huber_k08)

# --- Nearly OLS with large k ---
fit_huber_k5 <- fregre.huber(fd, y_contam, ncomp = 5, k = 5.0)
print(fit_huber_k5)
```

### The Tuning Parameter $k$

The parameter $k$ controls the transition between quadratic and linear
loss:

| $k$       | Behaviour        | Efficiency (Gaussian) | Robustness |
|-----------|------------------|-----------------------|------------|
| 0.5       | Very aggressive  | ~72%                  | High       |
| 1.0       | Moderate         | ~89%                  | Good       |
| **1.345** | **Default**      | **95%**               | **Good**   |
| 2.0       | Mild             | ~98%                  | Moderate   |
| $\infty$  | Identical to OLS | 100%                  | None       |

In practice, the default $k = 1.345$ is a good starting point. Decrease
$k$ if contamination is severe; increase $k$ if you trust most of the
data.

``` r
# Compare weights across k values
w_default <- fit_huber_contam$weights
w_aggressive <- fit_huber_k08$weights

cat("Default k=1.345:\n")
cat("  Mean weight (clean):", round(mean(w_default[-outlier_idx]), 3), "\n")
cat("  Mean weight (outlier):", round(mean(w_default[outlier_idx]), 3), "\n")

cat("Aggressive k=0.8:\n")
cat("  Mean weight (clean):", round(mean(w_aggressive[-outlier_idx]), 3), "\n")
cat("  Mean weight (outlier):", round(mean(w_aggressive[outlier_idx]), 3), "\n")
```

## Comparison: OLS vs L1 vs Huber

We compare all three methods on clean and contaminated data using a
train/test split.

``` r
# --- Train/test split ---
set.seed(123)
train_idx <- sample(n, 160)
test_idx <- setdiff(1:n, train_idx)

fd_train <- fd[train_idx, ]
fd_test <- fd[test_idx, ]

# ==========================================
# Clean data
# ==========================================
y_tr_clean <- y_clean[train_idx]
y_te_clean <- y_clean[test_idx]

fit_ols_c <- fregre.lm(fd_train, y_tr_clean, ncomp = 5)
fit_l1_c <- fregre.l1(fd_train, y_tr_clean, ncomp = 5)
fit_hub_c <- fregre.huber(fd_train, y_tr_clean, ncomp = 5)

pred_ols_c <- predict(fit_ols_c, fd_test)
pred_l1_c <- predict(fit_l1_c, fd_test)
pred_hub_c <- predict(fit_hub_c, fd_test)

# ==========================================
# Contaminated data
# ==========================================
y_tr_contam <- y_contam[train_idx]
y_te_contam <- y_contam[test_idx]

fit_ols_d <- fregre.lm(fd_train, y_tr_contam, ncomp = 5)
fit_l1_d <- fregre.l1(fd_train, y_tr_contam, ncomp = 5)
fit_hub_d <- fregre.huber(fd_train, y_tr_contam, ncomp = 5)

# Predict on test set --- compare to *clean* truth
pred_ols_d <- predict(fit_ols_d, fd_test)
pred_l1_d <- predict(fit_l1_d, fd_test)
pred_hub_d <- predict(fit_hub_d, fd_test)
```

### Results Table

``` r
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))
r2 <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)

results <- data.frame(
  Method = c("OLS (fregre.lm)", "L1 (fregre.l1)", "Huber (fregre.huber)"),
  RMSE_Clean = round(c(rmse(y_te_clean, pred_ols_c),
                        rmse(y_te_clean, pred_l1_c),
                        rmse(y_te_clean, pred_hub_c)), 4),
  R2_Clean = round(c(r2(y_te_clean, pred_ols_c),
                      r2(y_te_clean, pred_l1_c),
                      r2(y_te_clean, pred_hub_c)), 4),
  RMSE_Contam = round(c(rmse(y_te_clean, pred_ols_d),
                         rmse(y_te_clean, pred_l1_d),
                         rmse(y_te_clean, pred_hub_d)), 4),
  R2_Contam = round(c(r2(y_te_clean, pred_ols_d),
                       r2(y_te_clean, pred_l1_d),
                       r2(y_te_clean, pred_hub_d)), 4)
)

knitr::kable(results,
  col.names = c("Method", "RMSE", "R\u00b2", "RMSE", "R\u00b2"),
  caption = "Test-set performance on clean vs contaminated training data"
)
```

**Expected pattern:**

- On **clean data**, all three methods perform similarly. OLS is
  slightly more efficient because it is the maximum likelihood estimator
  under Gaussian errors.
- On **contaminated data**, OLS performance degrades sharply: the RMSE
  increases and $R^{2}$ drops. L1 and Huber remain close to their
  clean-data performance because they down-weight the outlier
  observations.

### Prediction with New Data

All three methods support the standard
[`predict()`](https://rdrr.io/r/stats/predict.html) generic:

``` r
# Predict from a robust model on new curves
new_curves <- fd[1:5, ]
pred_l1 <- predict(fit_l1_contam, new_curves)
pred_huber <- predict(fit_huber_contam, new_curves)

cat("L1 predictions:", round(pred_l1, 3), "\n")
cat("Huber predictions:", round(pred_huber, 3), "\n")
```

### Scalar Covariates

Both
[`fregre.l1()`](https://sipemu.github.io/fdars-r/reference/fregre.l1.md)
and
[`fregre.huber()`](https://sipemu.github.io/fdars-r/reference/fregre.huber.md)
support optional scalar covariates, just like
[`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md):

``` r
# Add a scalar covariate (e.g., sample temperature)
temperature <- rnorm(n, mean = 25, sd = 3)
y_with_temp <- y_clean + 0.5 * temperature

fit_l1_sc <- fregre.l1(fd, y_with_temp,
                        scalar.covariates = temperature, ncomp = 5)
print(fit_l1_sc)

# Prediction with scalar covariates
pred_sc <- predict(fit_l1_sc, fd[1:5, ],
                   scalar.covariates = temperature[1:5])
```

## When to Use Robust Regression

**Use
[`fregre.l1()`](https://sipemu.github.io/fdars-r/reference/fregre.l1.md)
when:**

- You suspect a substantial fraction (\> 5%) of response values may be
  corrupted by gross errors.
- You want maximum resistance to outliers and are willing to accept a
  small loss of efficiency on clean data.
- The error distribution is heavy-tailed (e.g., Cauchy or $t$ with few
  degrees of freedom).

**Use
[`fregre.huber()`](https://sipemu.github.io/fdars-r/reference/fregre.huber.md)
when:**

- You want a compromise: near-OLS efficiency on clean data with
  protection against moderate contamination.
- You are unsure whether outliers are present — Huber is a safe default
  that costs very little when the data are clean.
- You want to tune the robustness/efficiency trade-off via the parameter
  $k$.

**Stick with
[`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md)
when:**

- The data are clean and well-curated.
- You need the full downstream ecosystem: explainability (`explain()`),
  diagnostics (Cook’s distance, leverage, DFBETAS), uncertainty
  quantification (conformal prediction), and cross-validation
  ([`fregre.lm.cv()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.cv.md)).
- Efficiency is paramount and you trust the Gaussian error assumption.

### Limitations

- The robust methods protect against outliers in the **response** $y$.
  They do not address leverage points (outliers in the functional
  predictor $X(t)$). For functional outlier detection, see
  [`vignette("articles/outlier-detection")`](https://sipemu.github.io/fdars-r/articles/outlier-detection.md).
- The `fregre.robust` class does not currently support the full
  explainability and diagnostics ecosystem available for `fregre.lm`.
- Both robust methods rely on IRLS, which may not converge if the data
  are severely ill-conditioned. Check the `converged` field in the
  output.

## Mathematical Details

### FPC Projection

Both robust methods begin the same way as
[`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md):
perform FPCA on the functional predictor and project each curve onto $K$
principal components:

$$\xi_{ik} = \int\left( X_{i}(t) - \bar{X}(t) \right)\phi_{k}(t)\, dt,\quad k = 1,\ldots,K$$

This reduces the infinite-dimensional functional regression to a
finite-dimensional problem:

$$y_{i} = \alpha + \sum\limits_{k = 1}^{K}\gamma_{k}\xi_{ik} + \epsilon_{i}$$

The only difference from
[`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md)
is the loss function used to estimate
$\left( \alpha,\gamma_{1},\ldots,\gamma_{K} \right)$.

### IRLS Algorithm

Both L1 and Huber are solved by iteratively reweighted least squares. At
iteration $t$:

1.  Compute residuals:
    $r_{i}^{(t)} = y_{i} - \alpha^{(t)} - {\mathbf{ξ}}_{i}^{\top}{\mathbf{γ}}^{(t)}$
2.  Compute weights:
    $w_{i}^{(t)} = \psi\left( r_{i}^{(t)} \right)/r_{i}^{(t)}$
3.  Solve weighted least squares:
    $\left( \alpha^{(t + 1)},{\mathbf{γ}}^{(t + 1)} \right) = \arg\min\sum w_{i}^{(t)}\left( y_{i} - \alpha - {\mathbf{ξ}}_{i}^{\top}{\mathbf{γ}} \right)^{2}$

The weight function $\psi(r)/r$ differs by method:

- **L1:** $w_{i} = 1/\left| r_{i} \right|$ (with a floor to avoid
  division by zero)
- **Huber:** $w_{i} = \min\left( 1,k/\left| r_{i} \right| \right)$

Convergence is declared when the relative change in coefficients falls
below a tolerance, typically within 10–30 iterations.

## See Also

- [`vignette("articles/scalar-on-function")`](https://sipemu.github.io/fdars-r/articles/scalar-on-function.md)
  — overview of all scalar-on-function methods
- [`vignette("articles/outlier-detection")`](https://sipemu.github.io/fdars-r/articles/outlier-detection.md)
  — detecting outliers in functional data
- [`vignette("articles/regression-diagnostics")`](https://sipemu.github.io/fdars-r/articles/regression-diagnostics.md)
  — influence and diagnostics for
  [`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md)
- [`vignette("articles/explainability")`](https://sipemu.github.io/fdars-r/articles/explainability.md)
  — interpreting
  [`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md)
  models
- [`vignette("articles/cross-validation")`](https://sipemu.github.io/fdars-r/articles/cross-validation.md)
  — cross-validation framework

## References

- Huber, P.J. (1964). Robust Estimation of a Location Parameter. *Annals
  of Mathematical Statistics*, 35(1), 73–101.
- Huber, P.J. and Ronchetti, E.M. (2009). *Robust Statistics*, 2nd ed.
  Wiley.
- Maronna, R.A., Martin, R.D., Yohai, V.J., and Salibian-Barrera, M.
  (2019). *Robust Statistics: Theory and Methods*, 2nd ed. Wiley.
- Cardot, H., Crambes, C., and Sarda, P. (2005). Quantile Regression
  when the Covariates are Functions. *Journal of Nonparametric
  Statistics*, 17(7), 841–856.
