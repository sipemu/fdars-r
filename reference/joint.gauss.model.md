# Joint Amplitude-Phase Gaussian Generative Model

Fit a joint Gaussian generative model that captures both amplitude and
phase variability. Unlike
[`gauss.model`](https://sipemu.github.io/fdars-r/reference/gauss.model.md)
which models amplitude only, this model jointly models the amplitude and
warping components, allowing generated curves to exhibit realistic phase
variability.

## Usage

``` r
joint.gauss.model(
  karcher,
  n.components = 3,
  n.samples = 50,
  balance = 1,
  seed = 42
)
```

## Arguments

- karcher:

  An object of class 'karcher.mean' (result of
  [`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md)).

- n.components:

  Number of principal components to retain (default 3).

- n.samples:

  Number of random curves to generate (default 50).

- balance:

  Balance parameter controlling the relative weight of amplitude and
  phase components (default 1.0).

- seed:

  Random seed for reproducibility (default 42).

## Value

A list with components:

- samples:

  fdata of generated random curves

- eigenvalues:

  numeric vector of eigenvalues

- eigenfunctions:

  matrix of eigenfunctions

- scores:

  matrix of joint PCA scores

## References

Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
functional data using phase and amplitude separation. *Computational
Statistics & Data Analysis*, 61:50–66.

## See also

[`karcher.mean`](https://sipemu.github.io/fdars-r/reference/karcher.mean.md),
[`gauss.model`](https://sipemu.github.io/fdars-r/reference/gauss.model.md)

## Examples

``` r
# \donttest{
set.seed(1)
t <- seq(0, 1, length.out = 50)
X <- matrix(0, 20, 50)
for (i in 1:20) X[i, ] <- sin(2 * pi * t + runif(1, -0.3, 0.3))
fd <- fdata(X, argvals = t)
km <- karcher.mean(fd)
jgm <- joint.gauss.model(km, n.components = 2, n.samples = 10)
# }
```
