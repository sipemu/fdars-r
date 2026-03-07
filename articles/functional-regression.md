# Functional Regression: Scalar and Function Responses

## Introduction

Functional regression extends classical regression to settings where
predictors or responses are functions. **fdars** provides two
complementary directions:

- **Scalar-on-function regression**: scalar response, functional
  predictors
  ([`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md),
  [`functional.logistic()`](https://sipemu.github.io/fdars-r/reference/functional.logistic.md))
- **Function-on-scalar regression**: functional response, scalar
  predictors
  ([`fosr()`](https://sipemu.github.io/fdars-r/reference/fosr.md),
  [`fosr.fpc()`](https://sipemu.github.io/fdars-r/reference/fosr.fpc.md),
  [`fanova()`](https://sipemu.github.io/fdars-r/reference/fanova.md))

This article demonstrates both, including cross-validation and
hypothesis testing.

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

## Scalar-on-Function Regression

### FPC-Based Functional Linear Model

The functional linear model relates a scalar response $y_{i}$ to a
functional predictor $X_{i}(t)$ via:

$$y_{i} = \alpha + \int X_{i}(t)\beta(t)\, dt + \varepsilon_{i}$$

[`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md)
estimates $\beta(t)$ by projecting onto functional principal components.

``` r
set.seed(42)
n <- 80
m <- 100
t_grid <- seq(0, 1, length.out = m)

# Simulate functional predictor
X <- matrix(0, n, m)
for (i in 1:n) {
  X[i, ] <- sin(2 * pi * t_grid) * runif(1, 0.5, 2) +
    cos(4 * pi * t_grid) * rnorm(1, sd = 0.3) +
    rnorm(m, sd = 0.1)
}

# Scalar response depends on the first FPC direction
beta_true <- sin(2 * pi * t_grid)
y <- X %*% beta_true / m + rnorm(n, sd = 0.3)

fd <- fdata(X, argvals = t_grid)
fit <- fregre.lm(fd, y, ncomp = 3)
print(fit)
#> Functional Linear Model (FPC-based)
#> ====================================
#>   Number of observations: 80 
#>   Number of FPC components: 3 
#>   R-squared: 0.4371 
#>   Adjusted R-squared: 0.4149 
#>   GCV: 0.0968
```

### Cross-Validation for Component Selection

[`fregre.lm.cv()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.cv.md)
selects the optimal number of FPC components via k-fold
cross-validation.

``` r
cv_result <- fregre.lm.cv(fd, y, k.range = 1:8, nfold = 10)
cat("Optimal ncomp:", cv_result$optimal.k, "\n")
#> Optimal ncomp: 6
cat("CV errors:", round(cv_result$cv.errors, 4), "\n")
#> CV errors: 0.0957 0.0959 0.0971 0.0971 0.0978 0.0938 0.0942 0.0951
```

### Comparison with fregre.pc

[`fregre.lm()`](https://sipemu.github.io/fdars-r/reference/fregre.lm.md)
uses the Rust backend while
[`fregre.pc()`](https://sipemu.github.io/fdars-r/reference/fregre.pc.md)
uses the same FPC-based approach. Both produce similar fits:

``` r
fit_pc <- fregre.pc(fd, y, ncomp = 3)
cat("fregre.lm R-squared:", round(fit$r.squared, 4), "\n")
#> fregre.lm R-squared: 0.4371
cat("fregre.pc residual MSE:", round(mean(fit_pc$residuals^2), 4), "\n")
#> fregre.pc residual MSE: 0.0871
cat("fregre.lm residual MSE:", round(mean(fit$residuals^2), 4), "\n")
#> fregre.lm residual MSE: 0.0871
```

### Functional Logistic Regression

For binary classification,
[`functional.logistic()`](https://sipemu.github.io/fdars-r/reference/functional.logistic.md)
fits a logistic model with FPC features.

``` r
# Create binary response (0/1)
y_bin <- as.numeric(y > median(y))

logit_fit <- functional.logistic(fd, y_bin, ncomp = 3)
print(logit_fit)
#> Functional Logistic Regression
#> ==============================
#>   Number of observations: 80 
#>   Number of FPC components: 3 
#>   Accuracy: 0.85 
#>   Log-likelihood: -30.93 
#>   Iterations: 6

cat("Probabilities range:", range(logit_fit$probabilities), "\n")
#> Probabilities range: 0.01257656 0.9824597
```

## Function-on-Scalar Regression

### Penalized FOSR

When the response is a function and predictors are scalar, we have:

$$Y_{i}(t) = \mu(t) + \sum\limits_{j}x_{ij}\beta_{j}(t) + \varepsilon_{i}(t)$$

[`fosr()`](https://sipemu.github.io/fdars-r/reference/fosr.md) estimates
the coefficient functions $\beta_{j}(t)$ with optional ridge penalty.

``` r
set.seed(42)
n <- 60
m <- 80
t_grid <- seq(0, 1, length.out = m)

# Two scalar predictors: treatment and age
treatment <- rep(c(0, 1), each = n/2)
age <- rnorm(n)
predictors <- cbind(treatment, age)

# Functional response depends on predictors
Y <- matrix(0, n, m)
for (i in 1:n) {
  Y[i, ] <- sin(2 * pi * t_grid) +
    treatment[i] * 0.5 * cos(2 * pi * t_grid) +
    age[i] * 0.2 * t_grid +
    rnorm(m, sd = 0.15)
}

fd_y <- fdata(Y, argvals = t_grid)
fosr_fit <- fosr(fd_y, predictors, lambda = 0.1)
print(fosr_fit)
#> Function-on-Scalar Regression
#> =============================
#>   Number of observations: 60 
#>   Number of predictors: 2 
#>   Evaluation points: 80 
#>   R-squared: 0.6412 
#>   Lambda: 0.1
```

### Visualizing Coefficient Functions

``` r
plot(fosr_fit)
```

![](functional-regression_files/figure-html/fosr-plot-1.png)

### FPC-Based FOSR

[`fosr.fpc()`](https://sipemu.github.io/fdars-r/reference/fosr.fpc.md)
uses functional principal components instead of penalization:

``` r
fosr_fpc_fit <- fosr.fpc(fd_y, predictors, ncomp = 5)
cat("FPC-based R-squared:", round(fosr_fpc_fit$r.squared, 4), "\n")
#> FPC-based R-squared: 0.6421
```

### Prediction

``` r
new_pred <- matrix(c(1, 0.5,   # treated, average age
                     0, -1.0), # control, young
                   nrow = 2, byrow = TRUE)
pred <- predict(fosr_fit, new_pred)
```

## Functional ANOVA

[`fanova()`](https://sipemu.github.io/fdars-r/reference/fanova.md) tests
whether mean functions differ across groups using a permutation-based
F-test.

``` r
set.seed(42)
n_per_group <- 25
m <- 80
t_grid <- seq(0, 1, length.out = m)

# Three treatment groups with different mean functions
Y_anova <- matrix(0, 3 * n_per_group, m)
for (i in 1:n_per_group) {
  Y_anova[i, ] <- sin(2 * pi * t_grid) + rnorm(m, sd = 0.15)
  Y_anova[n_per_group + i, ] <- sin(2 * pi * t_grid) + 0.5 * cos(2 * pi * t_grid) +
    rnorm(m, sd = 0.15)
  Y_anova[2 * n_per_group + i, ] <- 2 * t_grid - 1 + rnorm(m, sd = 0.15)
}
groups <- rep(1:3, each = n_per_group)

fd_anova <- fdata(Y_anova, argvals = t_grid)
anova_result <- fanova(fd_anova, groups, n.perm = 500)
print(anova_result)
#> Functional ANOVA
#> ================
#>   Number of groups: 3 
#>   Number of observations: 75 
#>   Global F-statistic: 577.3765 
#>   P-value: 0.001996 
#>   Permutations: 500
```

### Visualizing Group Means

``` r
plot(anova_result)
```

![](functional-regression_files/figure-html/fanova-plot-1.png)

## See Also

- [`vignette("articles/regression")`](https://sipemu.github.io/fdars-r/articles/regression.md)
  for kernel and basis regression
- [`vignette("articles/functional-classification")`](https://sipemu.github.io/fdars-r/articles/functional-classification.md)
  for classification with
  [`fclassif()`](https://sipemu.github.io/fdars-r/reference/fclassif.md)
- [`vignette("articles/functional-mixed-models")`](https://sipemu.github.io/fdars-r/articles/functional-mixed-models.md)
  for repeated measures with
  [`fmm()`](https://sipemu.github.io/fdars-r/reference/fmm.md)
