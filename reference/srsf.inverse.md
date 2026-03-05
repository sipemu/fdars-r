# Inverse SRSF Transform

Reconstruct a curve from its SRSF representation: \\f(t) = f_0 +
\int_0^t q(s)\|q(s)\| ds\\.

## Usage

``` r
srsf.inverse(q, argvals, f0 = 0)
```

## Arguments

- q:

  Numeric vector of SRSF values.

- argvals:

  Numeric vector of evaluation points.

- f0:

  Initial value \\f_0\\ (default 0).

## Value

Numeric vector of reconstructed curve values.

## Examples

``` r
fd <- fdata(matrix(sin(seq(0, pi, length.out = 50)), 1, 50),
            argvals = seq(0, 1, length.out = 50))
q <- srsf.transform(fd)
f_hat <- srsf.inverse(q$data[1, ], fd$argvals, fd$data[1, 1])
```
