# Elastic Alignment for Functional Data

Functions for elastic curve alignment via SRSF transforms, Fisher-Rao
distance, and Karcher mean computation. SRSF Transform

## Usage

``` r
srsf.transform(fdataobj)
```

## Arguments

- fdataobj:

  An object of class 'fdata'.

## Value

An object of class 'fdata' containing the SRSF-transformed curves.

## Details

Compute the Square-Root Slope Function (SRSF) transformation of
functional data: \\q(t) = \text{sign}(f'(t)) \sqrt{\|f'(t)\|}\\.

## References

Srivastava, A., Klassen, E., Joshi, S.H., and Jermyn, I.H. (2011). Shape
analysis of elastic curves in Euclidean spaces. *IEEE Transactions on
Pattern Analysis and Machine Intelligence*, 33(7):1415–1428.

Tucker, J.D., Wu, W., and Srivastava, A. (2013). Generative models for
functional data using phase and amplitude separation. *Computational
Statistics & Data Analysis*, 61:50–66.

## Examples

``` r
fd <- fdata(matrix(rnorm(200), 20, 10), argvals = seq(0, 1, length.out = 10))
q <- srsf.transform(fd)
```
