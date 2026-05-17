# F7.1 — Thread-count policy for the OpenMP build.
#
# CRAN runs checks with a hard 2-thread budget, so the package must
# never spin up more threads than the user explicitly asks for. The
# resolved count is driven by the `mispitools.n_threads` option and
# defaults to 2 (CRAN-safe). `evaluate_evidence()` (F7.2) and the
# parallel kernels (F7.4) consume this helper; F7.1 only wires the
# policy and the OpenMP linkage probe.

#' Resolve the OpenMP thread count
#'
#' Returns the number of threads the engine should request, clamped to
#' the range supported by the build. Honours `getOption("mispitools.n_threads")`
#' and falls back to 2 to stay within the CRAN check budget. When the
#' package was built without OpenMP the result is always 1.
#'
#' @param requested Optional explicit thread count; overrides the option.
#' @return Integer, >= 1.
#' @noRd
mispi_n_threads <- function(requested = NULL) {
  info <- cpp_openmp_info(0L)
  if (!isTRUE(info$available)) {
    return(1L)
  }

  n <- requested
  if (is.null(n)) {
    n <- getOption("mispitools.n_threads", 2L)
  }
  n <- suppressWarnings(as.integer(n)[1L])
  if (is.na(n) || n < 1L) {
    n <- 1L
  }

  max_threads <- as.integer(info$max_threads)
  if (!is.na(max_threads) && max_threads >= 1L && n > max_threads) {
    n <- max_threads
  }
  n
}

#' Report the OpenMP build status
#'
#' @return List with `available` (logical), `max_threads`, and the
#'   `observed_threads` seen in a probe parallel region.
#' @noRd
mispi_openmp_status <- function() {
  cpp_openmp_info(mispi_n_threads())
}
