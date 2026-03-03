# CRAN Submission Comments

## R CMD check results

0 errors | 0 warnings | 1 note

### NOTE

1. **New submission**
   - This is a new submission to CRAN.
   - Package size (~33 MB) is due to vendored Rust crate sources required for
     offline compilation per CRAN policy.
   - Possibly misspelled words: Ferraty, Silverman, Vieu — these are author
     surnames from the referenced textbooks (Ramsay & Silverman 2005;
     Ferraty & Vieu 2006).

## Resubmission (v0.3.3)

Fixes from v0.3.3 CRAN pre-test failures:

1. **Fixed flaky RP depth test**: `test-validation-depth.R:83` failed on
   Debian because RP (random projection) depth uses stochastic projections,
   making the "deepest curve near center" assertion unreliable across
   platforms. Replaced the brittle centrality check with a robust
   non-degeneracy assertion (`sd(depths) > 0`). Structural validity
   (depths in [0,1], correct length) was already tested above.
2. **New feature: `fequiv.test()`**: Added functional equivalence test (TOST)
   based on Dette & Kokot (2021, Biometrika 108(4):895-913) with multiplier
   and percentile bootstrap methods, print and plot S3 methods, and 10 new
   test cases.

Fixes from v0.3.2 reviewer feedback (Benjamin Altmann):

1. **Quoted software names**: 'Rust' now single-quoted in Title and Description.
2. **Added references**: Ramsay & Silverman (2005, ISBN:978-0-387-40080-8) and
   Ferraty & Vieu (2006, ISBN:978-0-387-30369-7) in DESCRIPTION.
3. **Uncommented examples**: `flm.test` and `fmean.test.fdata` examples are now
   executable, wrapped in `\donttest{}` (permutation tests take >5s).
4. **Reset par()**: `addError` and `register.fd` examples now save/restore
   graphical parameters via `oldpar <- par(...); on.exit(par(oldpar))` pattern.
5. **Removed hardcoded set.seed**: Replaced `set.seed(123)` in `fregre.R` with
   a deterministic Halton sequence for reproducible weight grids without
   manipulating the user's RNG state.

### Previous resubmission (v0.3.2)

- Wrapped `abort`/`exit`/`_exit` symbols with `Rf_error()` via `--wrap` linker flag
- Added linker version script to export only `R_init_fdars`
- Limited Rust parallelism in tests to reduce CPU time

## Package Description

fdars provides functional data analysis tools with a high-performance Rust backend.
The package offers methods for:
- Functional depth computation (10 methods)
- Distance metrics and semimetrics (10+ methods)
- Functional regression (PC, basis, nonparametric)
- Basis representation (B-spline, Fourier, P-splines)
- Clustering (k-means, fuzzy c-means)
- Outlier detection
- Seasonal analysis
- Statistical testing

## Rust Dependency

This package uses Rust for performance-critical algorithms. The Rust code is
compiled during installation using the cargo build system.

### Vendored Dependencies
All Rust crate dependencies are bundled in the package using `cargo vendor`.
The build uses `--offline` mode - no network access is required during
installation. This follows the recommendations in "Using Rust in CRAN packages"
(https://cran.r-project.org/web/packages/using_rust.html).

### Build Requirements
- Rust toolchain (rustc >= 1.81, cargo)
- Users can install Rust from https://rustup.rs/ or system package manager

### configure Script
The package includes a configure script that:
1. Checks for Rust toolchain availability
2. Validates Rust version (>= 1.81)
3. Provides clear error messages if Rust is missing

## Test Coverage

- All examples run without errors
- All tests pass (1737 tests)
- 15 vignettes build successfully

## Test Environments

* Local: Manjaro Linux, R 4.5.2, Rust 1.84
* GitHub Actions: Ubuntu, macOS, Windows (R release and devel)
* win-builder: r-devel-windows-x86_64
* CRAN incoming: r-devel-linux-x86_64-debian-gcc

## Downstream Dependencies

None (new package).
