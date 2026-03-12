# CRAN Submission Comments

## R CMD check results

0 errors | 0 warnings | 1 note

### NOTE

1. **New submission**
   - This is a new submission to CRAN.
   - Possibly misspelled words: Ferraty, Silverman, Vieu — these are author
     surnames from the referenced textbooks (Ramsay & Silverman 2005;
     Ferraty & Vieu 2006).

## Changes since last submission (v0.3.3 → v0.5.0)

Major additions:
- Elastic alignment (Fisher-Rao framework): SRSF, dynamic programming alignment,
  Karcher mean, amplitude/phase decomposition
- Functional regression: function-on-scalar, functional ANOVA, functional
  mixed models, functional logistic regression
- Functional classification: LDA, QDA, kNN, kernel, DD-plot with cross-validation
- GMM clustering with automatic K selection
- Andrews transformation for multivariate-to-functional projection
- Tolerance bands (FPCA, conformal, SCB, exponential, elastic)
- Landmark registration and TSRVF
- Simulation toolbox (Karhunen-Loève, irregular/sparse data)
- Streaming depth computation
- Soft-DTW distance and barycenter

CRAN compliance improvements:
- Rayon thread count limited to 2 via RAYON_NUM_THREADS in .onLoad()
- inst/COPYRIGHTS file with all vendored Rust crate authors and licenses
- Copyright field added to DESCRIPTION
- SystemRequirements specifies rustc version (>= 1.81)
- Portable Makefiles (removed GNU extensions)
- Removed hidden files and long paths from vendored crates
- Internal Rust wrapper functions no longer generate man pages

## Rust Dependency

This package uses Rust for performance-critical algorithms. The Rust code is
compiled during installation using the cargo build system.

### Vendored Dependencies
All Rust crate dependencies are bundled in the package using `cargo vendor`.
The build uses `--offline` mode — no network access is required during
installation. This follows the recommendations in "Using Rust in CRAN packages"
(https://cran.r-project.org/web/packages/using_rust.html).

### Thread Limiting
The package sets `RAYON_NUM_THREADS=2` in `.onLoad()` to comply with CRAN's
policy on parallel resource usage. Users can override by setting the environment
variable before loading the package.

### Symbol Visibility
On Linux, `--wrap=abort,exit,_exit` linker flags and a `--version-script`
hide Rust stdlib symbols (abort, exit, _exit) that R CMD check would flag.
On macOS, `-exported_symbols_list` achieves the same by exporting only
`R_init_fdars`. These symbols are never actually called by the package.

### Panic Safety
All Rust functions called from R are wrapped by extendr's `catch_unwind`
mechanism, which converts Rust panics into R errors via `Rf_error()`.

### Build Requirements
- Rust toolchain (rustc >= 1.81, cargo)
- Users can install Rust from https://rustup.rs/ or system package manager

### configure Script
The package includes a configure script that:
1. Checks for Rust toolchain availability with helpful install instructions
2. Validates Rust version (>= 1.81)
3. Restores vendored cargo checksum files
4. Generates platform-appropriate Makevars

## Test Coverage

- All examples run without errors
- 842 test cases across 31 test files
- 31 vignette articles build successfully

## Test Environments

* Local: Manjaro Linux, R 4.5.2, Rust 1.93
* GitHub Actions: Ubuntu (R release + devel), macOS (R release), Windows (R release)

## Downstream Dependencies

None (new package).
