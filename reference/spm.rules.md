# Evaluate Control Chart Rules

Applies Western Electric or Nelson run rules to a sequence of control
chart values. Identifies patterns that indicate non-random behavior
(e.g., runs, trends, points beyond control limits).

## Usage

``` r
spm.rules(values, center, sigma, rule.set = c("western.electric", "nelson"))
```

## Arguments

- values:

  Numeric vector of chart values (in time order).

- center:

  Center line value (e.g., process mean).

- sigma:

  Standard deviation (for zone boundaries at 1, 2, 3 sigma).

- rule.set:

  Which rule set to apply: `"western.electric"` (4 rules) or `"nelson"`
  (8 rules).

## Value

An object of class `spm.rules` with components:

- rules:

  Character vector of violated rule names

- indices:

  List of integer vectors; each gives the observation indices
  (1-indexed) involved in that violation

- values:

  The input values

- center:

  The center line

- sigma:

  The standard deviation

- rule.set:

  The rule set used

## Examples

``` r
# \donttest{
set.seed(1)
values <- c(rnorm(20), rnorm(5, mean = 3))
result <- spm.rules(values, center = 0, sigma = 1,
                    rule.set = "western.electric")
result
#> SPM Control Chart Rule Evaluation
#>   Rule set: western.electric 
#>   Observations: 25 
#>   Center: 0 
#>   Sigma: 1 
#>   Violations found: 11 
#>     WE1: indices [21]
#>     WE1: indices [22]
#>     WE1: indices [23]
#>     WE1: indices [25]
#>     WE2: indices [21, 22]
#>     WE2: indices [21, 22, 23]
#>     WE2: indices [22, 23]
#>     WE2: indices [23, 25]
#>     WE3: indices [21, 22, 23, 24]
#>     WE3: indices [21, 22, 23, 24, 25]
#>     WE4: indices [18, 19, 20, 21, 22, 23, 24, 25]
# }
```
