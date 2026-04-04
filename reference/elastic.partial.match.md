# Elastic Partial Curve Matching

Find the optimal partial match of a template curve within a longer
target curve using the elastic framework. This identifies the
sub-interval of the target that best matches the template, along with
the optimal warping.

## Usage

``` r
elastic.partial.match(
  template,
  target,
  argvals.template,
  argvals.target,
  lambda = 0,
  min.span = 0.3
)
```

## Arguments

- template:

  Numeric vector of the template curve's values.

- target:

  Numeric vector of the target curve's values.

- argvals.template:

  Numeric vector of evaluation points for the template.

- argvals.target:

  Numeric vector of evaluation points for the target.

- lambda:

  Regularisation parameter controlling warping smoothness (default 0).

- min.span:

  Minimum fraction of the target domain that the matched sub-interval
  must span (default 0.3).

## Value

A list with components:

- gamma:

  numeric vector of the warping function on the matched sub-interval

- match_start:

  start of the matched sub-interval on the target domain

- match_end:

  end of the matched sub-interval on the target domain

- distance:

  elastic distance of the partial match

- target_aligned:

  numeric vector of the matched and aligned target segment

## References

Srivastava, A. and Klassen, E. (2016). *Functional and Shape Data
Analysis*. Springer.

## Examples

``` r
# \donttest{
t1 <- seq(0, 1, length.out = 50)
t2 <- seq(0, 2, length.out = 100)
template <- sin(2 * pi * t1)
target <- c(rep(0, 20), sin(2 * pi * seq(0, 1, length.out = 50)), rep(0, 30))
res <- elastic.partial.match(template, target, t1, t2)
# }
```
