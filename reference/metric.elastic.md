# Elastic Distance (Metric Dispatcher Alias)

Alias for
[`elastic.distance`](https://sipemu.github.io/fdars-r/reference/elastic.distance.md)
for use with the
[`metric`](https://sipemu.github.io/fdars-r/reference/metric.md)
dispatcher.

## Usage

``` r
metric.elastic(fdataobj, fdataref = NULL, ...)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

- fdataref:

  Optional reference 'fdata'. If NULL, computes self-distances.

- ...:

  Additional arguments (ignored).

## Value

A distance matrix (numeric matrix).

## Examples

``` r
# \donttest{
fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
D <- metric.elastic(fd)
# }
```
