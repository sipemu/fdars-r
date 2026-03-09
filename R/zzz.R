#' Package initialization
#'
#' @useDynLib fdars, .registration = TRUE
#' @importFrom methods is
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # CRAN policy: never use more than 2 simultaneous threads/cores
  # Rayon (Rust parallelism library) defaults to all available CPUs
  if (Sys.getenv("RAYON_NUM_THREADS") == "") {
    Sys.setenv(RAYON_NUM_THREADS = "2")
  }
}

.onUnload <- function(libpath) {
  library.dynam.unload("fdars", libpath)
}
