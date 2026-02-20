# CRAN Submission Comments

## R CMD check results

0 errors | 0 warnings | 1 note

### NOTE

1. **New submission**
   - This is a new submission to CRAN.
   - Package size (~33 MB) is due to vendored Rust crate sources required for
     offline compilation per CRAN policy.

## Resubmission

Fixes from v0.3.1 submission feedback:

- **Compiled code WARNING (exit/abort)**: Used linker `--wrap` flag to intercept
  `abort`/`exit`/`_exit` from the Rust standard library, replacing them with
  `Rf_error()` wrappers that safely convert to R errors. Combined with a linker
  version script that hides all symbols except `R_init_fdars`.
- **Test CPU time NOTE**: Limited Rust thread pool to 2 threads via
  `RAYON_NUM_THREADS=2` in test setup.

## Changes Since Last Submission

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
- All tests pass (1714 tests)
- 15 vignettes build successfully

## Test Environments

* Local: Manjaro Linux, R 4.5.2, Rust 1.84
* GitHub Actions: Ubuntu, macOS, Windows (R release and devel)

## Downstream Dependencies

None (new package).
