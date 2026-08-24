## Submission

This is a major release (1.4.0 -> 2.0.0). It adds an exact C++ engine for
likelihood-ratio distributions alongside the existing simulation workflow.
All functions from the 1.x series remain exported and unchanged in behaviour.

## Test environments

* local: Ubuntu 24.04, R 4.5.2, gcc 13 (`R CMD check --as-cran`)

## R CMD check results

0 errors | 0 warnings | 0 notes attributable to the package.

The local run reports one warning and two notes, all three of which come from
the checking environment rather than from the package:

* `'qpdf' is needed for checks on size reduction of PDFs` — qpdf is not
  installed on the machine used for the check.
* `Skipping checking HTML validation: no command 'tidy' found` — likewise.
* `Compilation used the following non-portable flag(s): '-mno-omit-leaf-frame-pointer'`
  — this flag is injected by the Debian/Ubuntu build of R through
  `/usr/lib/R/etc/Makeconf`; it does not appear in the package. `src/Makevars`
  sets only `CXX_STD = CXX17`, `-I.`, and `$(SHLIB_OPENMP_CXXFLAGS)`.

## Compiled code

The package contains a C++17 kernel reached through Rcpp and RcppArmadillo.
OpenMP is requested through `$(SHLIB_OPENMP_CXXFLAGS)` so that R resolves the
correct flag for the toolchain, and every region touching the runtime is
guarded by `#ifdef _OPENMP`. The package builds and runs single-threaded when
OpenMP is unavailable, with identical results.

## Package size

The previous development tarball carried a 6.6 MB tutorial video under
`man/figures/`. It has been excluded from the build; the tutorial is linked
from the README to its hosted location. The source tarball is now 1.7 MB.

## Downstream dependencies

There are no reverse dependencies on CRAN.
