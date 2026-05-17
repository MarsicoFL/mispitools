## F7.1 — OpenMP portability + linkage.
##
## The package must build and run whether or not the active toolchain
## ships an OpenMP runtime (src/Makevars uses $(SHLIB_OPENMP_CXXFLAGS),
## which is empty when OpenMP is absent, and every omp_* call is guarded
## by #ifdef _OPENMP). These tests pin the contract of the linkage probe
## and the CRAN-safe thread-count policy without assuming OpenMP is
## present, while still *proving* the runtime is linked (not just that
## the macro is defined) when it is.
##
## Options are saved/restored with base R (no withr) to avoid an
## undeclared test dependency under R CMD check.

with_n_threads_option <- function(value, code) {
  old <- options(mispitools.n_threads = value)
  on.exit(options(old), add = TRUE)
  force(code)
}

test_that("cpp_openmp_info has a stable contract", {
  info <- mispitools:::cpp_openmp_info(0L)
  expect_named(info, c("available", "max_threads", "observed_threads"))
  expect_type(info$available, "logical")
  expect_length(info$available, 1L)
  expect_gte(as.integer(info$max_threads), 1L)
  expect_gte(as.integer(info$observed_threads), 1L)
})

test_that("OpenMP runtime is actually linked when available", {
  info <- mispitools:::cpp_openmp_info(0L)
  if (!isTRUE(info$available)) {
    skip("package built without OpenMP (SHLIB_OPENMP_CXXFLAGS empty)")
  }
  ## A real parallel region ran: the omp runtime is linked, not just
  ## the _OPENMP macro defined at compile time.
  expect_gte(as.integer(info$max_threads), 1L)
  expect_gte(as.integer(info$observed_threads), 1L)

  ## num_threads clause must cap the region. With >=2 threads available
  ## a request for 2 must be honoured exactly.
  if (as.integer(info$max_threads) >= 2L) {
    capped <- mispitools:::cpp_openmp_info(2L)
    expect_equal(as.integer(capped$observed_threads), 2L)
  }
  ## Requesting a single thread must always serialise the region.
  one <- mispitools:::cpp_openmp_info(1L)
  expect_equal(as.integer(one$observed_threads), 1L)
})

test_that("mispi_n_threads defaults to the CRAN-safe budget", {
  info <- mispitools:::cpp_openmp_info(0L)
  with_n_threads_option(NULL, {
    n <- mispitools:::mispi_n_threads()
    if (isTRUE(info$available)) {
      expect_equal(n, min(2L, as.integer(info$max_threads)))
    } else {
      expect_equal(n, 1L)
    }
  })
})

test_that("mispi_n_threads honours the option and clamps to the build", {
  info <- mispitools:::cpp_openmp_info(0L)
  max_threads <- as.integer(info$max_threads)

  with_n_threads_option(1L, {
    expect_equal(mispitools:::mispi_n_threads(), 1L)
  })

  ## Absurd request is clamped to what the build supports (1 when no
  ## OpenMP, max_threads otherwise).
  with_n_threads_option(9999L, {
    n <- mispitools:::mispi_n_threads()
    if (isTRUE(info$available)) {
      expect_equal(n, max_threads)
    } else {
      expect_equal(n, 1L)
    }
  })

  ## Garbage / non-positive falls back to a single thread.
  with_n_threads_option(-3L, {
    expect_equal(mispitools:::mispi_n_threads(), 1L)
  })
  expect_equal(mispitools:::mispi_n_threads(requested = 1L), 1L)
})

test_that("no-OpenMP fallback degrades to single-thread cleanly", {
  info <- mispitools:::cpp_openmp_info(0L)
  if (isTRUE(info$available)) {
    skip("OpenMP available; fallback branch exercised on no-OpenMP builds")
  }
  expect_false(info$available)
  expect_equal(as.integer(info$max_threads), 1L)
  expect_equal(as.integer(info$observed_threads), 1L)
  expect_equal(mispitools:::mispi_n_threads(), 1L)
})
