# Functional Mixed Models for Repeated Measures

## Introduction

In many studies, each subject is observed multiple times, producing
repeated functional measurements. A **functional mixed model** (FMM)
decomposes the variation into:

$$Y_{ij}(t) = \mu(t) + \sum\limits_{k}x_{ij,k}\beta_{k}(t) + b_{i}(t) + \varepsilon_{ij}(t)$$

where:

- $\mu(t)$ is the population mean function
- $\beta_{k}(t)$ are fixed-effect coefficient functions (e.g., treatment
  effects)
- $b_{i}(t)$ are subject-level random effects
- $\varepsilon_{ij}(t)$ is measurement error

**fdars** provides
[`fmm()`](https://sipemu.github.io/fdars-r/reference/fmm.md) for
fitting,
[`fmm.predict()`](https://sipemu.github.io/fdars-r/reference/fmm.predict.md)
for prediction, and
[`fmm.test.fixed()`](https://sipemu.github.io/fdars-r/reference/fmm.test.fixed.md)
for testing fixed effects via permutation.

``` r
library(fdars)
#> 
#> Attaching package: 'fdars'
#> The following objects are masked from 'package:stats':
#> 
#>     cov, decompose, deriv, median, sd, var
#> The following object is masked from 'package:base':
#> 
#>     norm
library(ggplot2)
theme_set(theme_minimal())
```

## Why Mixed Models for Functional Data?

When subjects contribute multiple curves, the observations are no longer
independent. Ignoring this structure has serious consequences:

1.  **Pseudoreplication.** Treating all $n_{\text{total}}$ curves as
    independent inflates the effective sample size. Standard errors are
    too small, p-values are too optimistic, and Type I error rates
    exceed their nominal level.

2.  **Confounding sources of variation.** Without subject-level random
    effects, between-subject variability is absorbed into the residual,
    blurring the distinction between “subjects differ” and “measurements
    are noisy.”

3.  **Inefficient estimation.** The within-subject correlation carries
    information about fixed effects. Mixed models exploit this
    correlation structure, producing more efficient (lower-variance)
    estimates than methods that ignore it.

The functional mixed model resolves these issues by decomposing each
curve into four additive components:

| Component       | Notation                       | Interpretation                                       |
|-----------------|--------------------------------|------------------------------------------------------|
| Population mean | $\mu(t)$                       | Average curve across all subjects and conditions     |
| Fixed effects   | $\sum_{k}x_{ij,k}\beta_{k}(t)$ | Systematic effects of covariates (e.g., treatment)   |
| Random effects  | $b_{i}(t)$                     | Subject-specific deviations from the population mean |
| Residual        | $\varepsilon_{ij}(t)$          | Within-subject measurement noise                     |

This decomposition is the functional analog of the classical scalar
linear mixed model
$y_{ij} = \mu + x_{ij}^{T}\beta + b_{i} + \varepsilon_{ij}$, with each
scalar parameter replaced by a function of $t$.

## Mathematical Framework

### Model Specification

The full model for the $j$-th observation from subject $i$ is:

$$Y_{ij}(t) = \mu(t) + \sum\limits_{p = 1}^{P}\beta_{p}(t)\, x_{ip} + b_{i}(t) + \varepsilon_{ij}(t)$$

with distributional assumptions:

- $b_{i}(t) \sim \text{GP}\left( 0,\sigma_{u}^{2} \right)$:
  subject-specific random effects, independent across subjects
- $\varepsilon_{ij}(t) \sim \text{GP}\left( 0,\sigma_{\varepsilon}^{2} \right)$:
  residual noise, independent across observations

Here GP denotes a Gaussian process. The key parameters are the two
variance components $\sigma_{u}^{2}$ (between-subject) and
$\sigma_{\varepsilon}^{2}$ (within-subject).

### FPC-Based Estimation

Direct estimation of the functional parameters
$\mu(t),\beta_{p}(t),b_{i}(t)$ is infinite-dimensional. The FPC approach
reduces this to a set of scalar mixed models.

**Step 1.** Compute the FPC decomposition of the pooled data:

$$Y_{ij}(t) - \bar{Y}(t) \approx \sum\limits_{k = 1}^{K}\xi_{ij,k}\,\phi_{k}(t)$$

where $\phi_{k}(t)$ are the eigenfunctions and $\xi_{ij,k}$ are the
scores.

**Step 2.** For each component $k = 1,\ldots,K$, fit a scalar mixed
model:

$$\xi_{ij,k} = \gamma_{0k} + \sum\limits_{p = 1}^{P}\gamma_{pk}\, x_{ip} + u_{ik} + e_{ijk}$$

where $u_{ik} \sim \mathcal{N}\left( 0,\sigma_{u,k}^{2} \right)$ are
subject-level random effects and
$e_{ijk} \sim \mathcal{N}\left( 0,\sigma_{\varepsilon,k}^{2} \right)$
are residuals.

**Step 3.** Reconstruct the functional parameters:

$$\widehat{\mu}(t) = \bar{Y}(t) + \sum\limits_{k}{\widehat{\gamma}}_{0k}\,\phi_{k}(t),\qquad{\widehat{\beta}}_{p}(t) = \sum\limits_{k}{\widehat{\gamma}}_{pk}\,\phi_{k}(t)$$

$${\widehat{b}}_{i}(t) = \sum\limits_{k}{\widehat{u}}_{ik}\,\phi_{k}(t)$$

This approach is elegant: the infinite-dimensional functional mixed
model decomposes into $K$ independent scalar mixed models, each of which
can be fitted using standard REML estimation.

### REML Estimation

Within each FPC component, the variance components are estimated via
**Restricted Maximum Likelihood (REML)**, which accounts for the loss of
degrees of freedom from estimating fixed effects.

The key quantities are:

**Shrinkage weights:**

$$w_{i} = \frac{\sigma_{u}^{2}}{\sigma_{u}^{2} + \sigma_{\varepsilon}^{2}/n_{i}}$$

where $n_{i}$ is the number of observations for subject $i$. When
$n_{i}$ is large or $\sigma_{u}^{2} \gg \sigma_{\varepsilon}^{2}$,
$w_{i} \approx 1$ and the random effect estimate is close to the raw
subject mean. When $n_{i}$ is small or
$\sigma_{u}^{2} \ll \sigma_{\varepsilon}^{2}$, $w_{i} \approx 0$ and the
estimate is shrunk toward zero (the population mean).

**BLUP (Best Linear Unbiased Predictor) for random effects:**

$${\widehat{u}}_{i} = w_{i}({\bar{\xi}}_{i} - {\bar{x}}_{i}^{T}\widehat{\gamma})$$

The BLUP “borrows strength” across subjects: subjects with few
observations are pulled toward the population mean, while subjects with
many observations retain their individual character. This is the
functional analog of the well-known shrinkage property of scalar mixed
models.

**REML divisor:** The variance estimate uses $n - p$ in the denominator
(where $p$ is the number of fixed-effect parameters) rather than $n$,
correcting for the bias that would otherwise arise from estimating fixed
effects and variance components simultaneously.

### Hypothesis Testing for Fixed Effects

Testing whether a covariate has a significant functional effect:

$$H_{0}:\beta_{p}(t) = 0{\mspace{6mu}\text{for all}\mspace{6mu}}t$$

is analogous to functional ANOVA. The test uses a permutation approach
that respects the within-subject correlation:

1.  Compute the observed F-statistic $F_{\text{obs}}$ from the fitted
    model.
2.  For $b = 1,\ldots,B$: permute subject-level labels (not individual
    observations) and refit the model to obtain $F_{b}^{*}$.
3.  The p-value is:

$$p = \frac{1 + \#\{ F_{b}^{*} \geq F_{\text{obs}}\}}{B + 1}$$

Permuting at the subject level (rather than the observation level)
preserves the within-subject correlation structure under the null
hypothesis, ensuring valid inference even with complex dependency
patterns.

## Simulated Repeated Measures Data

We simulate a clinical trial with 15 subjects, each measured at 3 time
points, with a binary treatment covariate.

``` r
set.seed(42)
n_subjects <- 15
n_visits <- 3
n_total <- n_subjects * n_visits
m <- 80
t_grid <- seq(0, 1, length.out = m)

# Subject-specific random effects
subject_effects <- matrix(0, n_subjects, m)
for (i in 1:n_subjects) {
  subject_effects[i, ] <- 0.3 * rnorm(1) * sin(2 * pi * t_grid) +
    0.2 * rnorm(1) * cos(4 * pi * t_grid)
}

# Treatment assignment (first 8 subjects treated, rest control)
treatment <- c(rep(1, 8), rep(0, 7))

# Generate observations
Y <- matrix(0, n_total, m)
subject_ids <- rep(1:n_subjects, each = n_visits)
covariates <- matrix(0, n_total, 1)

for (i in 1:n_total) {
  subj <- subject_ids[i]
  trt <- treatment[subj]
  covariates[i, 1] <- trt

  Y[i, ] <- sin(2 * pi * t_grid) +                         # Mean function
    trt * 0.5 * cos(2 * pi * t_grid) +                      # Treatment effect
    subject_effects[subj, ] +                                 # Random effect
    rnorm(m, sd = 0.1)                                        # Noise
}

fd <- fdata(Y, argvals = t_grid)
```

The simulation has three variance sources: the treatment effect
($0.5\cos(2\pi t)$ for treated subjects), subject-level random effects
(smooth random curves with $\sigma \approx 0.3$), and measurement noise
($\sigma = 0.1$). The model should recover these components.

## Fitting the Functional Mixed Model

``` r
fit <- fmm(fd, subject_ids, covariates = covariates, ncomp = 3)
print(fit)
#> Functional Mixed Model
#> ======================
#>   Number of observations: 45 
#>   Number of subjects: 15 
#>   FPC components: 3 
#>   Residual variance (sigma2_eps): 0.000159 
#>   Fixed effect covariates: 1
```

The output shows the estimated variance components and the number of FPC
components used. The ratio $\sigma_{u}^{2}/\sigma_{\varepsilon}^{2}$
(the intraclass correlation in functional form) indicates how much of
the residual variation is between-subject vs. within-subject.

## Visualizing Fixed and Random Effects

The plot shows the population mean (red) overlaid with subject-specific
random effects (blue):

``` r
plot(fit)
```

![](functional-mixed-models_files/figure-html/fmm-plot-1.png)

Subject-specific curves that deviate substantially from the population
mean indicate large random effects for those subjects. The spread of the
blue curves reflects $\sigma_{u}^{2}$: wider spread means more
between-subject variability.

## Examining Model Components

``` r
# Residual variance
cat("Residual variance:", round(fit$sigma2.eps, 4), "\n")
#> Residual variance: 2e-04

# Number of subjects
cat("Number of subjects:", fit$n.subjects, "\n")
#> Number of subjects: 15

# FPC components used for random effects
cat("FPC components:", fit$ncomp, "\n")
#> FPC components: 3
```

The residual variance ${\widehat{\sigma}}_{\varepsilon}^{2}$ should be
close to $0.1^{2} = 0.01$ (the simulation noise level). A much larger
value would suggest the model is not capturing all systematic variation.

## Population-Level Prediction

[`fmm.predict()`](https://sipemu.github.io/fdars-r/reference/fmm.predict.md)
generates predictions using fixed effects only (no random effects),
useful for predicting the population-level response for new covariate
values.

``` r
# Predict for treated and control groups
new_cov <- matrix(c(1, 0), nrow = 2)
pred <- fmm.predict(fit, new.covariates = new_cov)

# Plot predictions
df_pred <- data.frame(
  t = rep(t_grid, 2),
  value = as.vector(t(pred$data)),
  group = factor(rep(c("Treated", "Control"), each = m))
)

ggplot(df_pred, aes(x = t, y = value, color = group)) +
  geom_line(linewidth = 1) +
  labs(x = "t", y = "Y(t)",
       title = "Predicted Population-Level Responses",
       color = "Group")
```

![](functional-mixed-models_files/figure-html/predict-1.png)

The gap between the treated and control curves is the estimated
treatment effect ${\widehat{\beta}}_{1}(t)$. Where the curves diverge
most is where the treatment has its strongest functional effect.

## Testing Fixed Effects

[`fmm.test.fixed()`](https://sipemu.github.io/fdars-r/reference/fmm.test.fixed.md)
tests whether the treatment has a significant effect using a permutation
F-test that respects within-subject correlation:

``` r
test_result <- fmm.test.fixed(fd, subject_ids, covariates,
                               ncomp = 3, n.perm = 500, seed = 123)
print(test_result)
#> Fixed Effects Permutation Test
#> ==============================
#>   Permutations: 500 
#> 
#>             F.statistic  P.value
#> Covariate 1      0.0394 0.001996
```

A small p-value indicates that the treatment effect $\beta_{1}(t)$ is
significantly different from zero. Because the permutation is done at
the subject level, this test correctly accounts for the repeated
measures structure — unlike a naive test that would treat all 45 curves
as independent.

## Comparison: Ignoring Subject Structure

What happens if we ignore the repeated measures structure and treat all
observations as independent?

``` r
# Naive approach: function-on-scalar regression ignoring subjects
fosr_naive <- fosr(fd, covariates)
cat("Naive FOSR R-squared:", round(fosr_naive$r.squared, 4), "\n")
#> Naive FOSR R-squared: 0.3143

# FMM accounts for subject-level variation
cat("FMM residual variance:", round(fit$sigma2.eps, 4), "\n")
#> FMM residual variance: 2e-04
```

The FMM separates subject-level variation from residual noise, giving a
cleaner estimate of the treatment effect. The naive FOSR lumps both
sources together, producing an inflated residual and potentially
misleading inference (standard errors that are too small because the
effective sample size is overstated).

## Interpretation and Diagnostics

### Variance Ratio

The ratio
${\widehat{\sigma}}_{u}^{2}/\left( {\widehat{\sigma}}_{u}^{2} + {\widehat{\sigma}}_{\varepsilon}^{2} \right)$
is the **intraclass correlation coefficient (ICC)**. It describes the
proportion of non-fixed-effect variability attributable to subject
differences:

- **ICC near 1**: subjects are very different from each other;
  within-subject replicates are nearly identical. Mixed models are
  essential.
- **ICC near 0**: subjects are similar; most variation is within-subject
  noise. A simpler model (ignoring subjects) may suffice.

### Shrinkage Toward the Population Mean

Random effects are shrunk toward zero by the factor $w_{i}$. Subjects
with fewer observations are shrunk more heavily. This is desirable — it
prevents overfitting to individual subjects with sparse data — but can
be misleading if interpreted as “this subject truly has a small effect.”
The amount of shrinkage is itself informative: heavy shrinkage indicates
that the data for that subject are too noisy to estimate a reliable
individual effect.

### When to Suspect Model Misspecification

- **Systematic residual patterns**: if residuals show temporal structure
  (e.g., autocorrelation), the residual process $\varepsilon_{ij}(t)$
  may not be white noise. Consider adding more FPC components or
  modeling serial correlation.
- **Non-constant variance**: if residual variance changes across the
  domain, a heteroscedastic model may be needed.
- **Random effects not Gaussian**: the REML estimator is still
  consistent under non-Gaussianity, but prediction intervals for random
  effects may be poorly calibrated.

## References

- Morris, J.S. (2015). Functional regression. *Annual Review of
  Statistics and Its Application*, 2, 321–359.
- Di, C.Z., Crainiceanu, C.M., Caffo, B.S., and Punjabi, N.M. (2009).
  Multilevel functional principal component analysis. *Annals of Applied
  Statistics*, 3(1), 458–488.
- Greven, S., Crainiceanu, C., Caffo, B., and Reich, D. (2010).
  Longitudinal functional principal component analysis. *Electronic
  Journal of Statistics*, 4, 1022–1054.
- Crainiceanu, C.M., Staicu, A.M., and Di, C.Z. (2009). Generalized
  multilevel functional regression. *Journal of the American Statistical
  Association*, 104(488), 1550–1561.

## See Also

- [`vignette("articles/functional-regression")`](https://sipemu.github.io/fdars-r/articles/functional-regression.md)
  for scalar-on-function and function-on-scalar regression
- [`vignette("articles/fpca")`](https://sipemu.github.io/fdars-r/articles/fpca.md)
  for the FPC decomposition underlying the random effects
