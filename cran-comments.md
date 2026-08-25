## Resubmission

This is a resubmission of version 2.0.0. The 2026-08-24 submission was
rejected by the incoming auto-check for

    Overall checktime 17 min > 10 min

on r-devel-windows-x86_64, of which the test suite accounted for 12 minutes.

The fix is confined to the test suite. Relative to the rejected tarball,
`R/`, `src/`, `man/*.Rd`, `NAMESPACE` and `DESCRIPTION` are unchanged, so
the exported API and the compiled kernel are identical. The only other
change is `README.md`, which gained a section reporting measured
comparisons of the exact engine against the simulation workflow, with
four figures under `man/figures/`.

### What was slow, and what was done about it

The cost sat in the cross-engine verification files (`tests/testthat/
test-vs-*.R`). Those compare the package engines against independent
oracles -- pedprobr, pedmut, Familias and forrel -- over a full grid of
pedigree x mutation model x marker, plus two Monte Carlo checks at
N = 20000 simulated profiles. The slowest cells are dominated by the
oracle's brute-force enumeration over 5- and 8-member pedigrees, not by
mispitools itself.

Following the suggestion in the rejection message, each grid is now split
in two: a representative subset that always runs, and the exhaustive
remainder gated on an environment variable. The gate is
`testthat::skip_on_cran()`, i.e. the remainder runs only when `NOT_CRAN`
is set to "true", which is the case under `devtools::test()`,
`devtools::check()` and CI, and is not the case on CRAN. The rationale
and the split are documented in `tests/testthat/helper-cran.R` and at
each gate.

Concretely:

* Toy data throughout: the always-on cells keep the small pedigrees and
  the smaller allele sets (3 alleles rather than 4, 2 rather than 3),
  which is where the state space of the comparison actually lives.
* Fewer iterations: the two forrel Monte Carlo checks run at N = 2000 on
  CRAN instead of N = 20000. Their gates are stated as multiples of the
  sample standard error of the mean, which is recomputed from the draws,
  so they stay correctly calibrated at the smaller N.
* Conditional tests: the redundant cells of each grid -- larger allele
  sets, additional pedigree topologies, interior recombination fractions,
  identities already covered by another always-on cell -- are skipped
  unless `NOT_CRAN` is set.

No test was deleted, and every code path that had an oracle check still
has one on CRAN.

### Measured effect

`R CMD check --as-cran` on the maintainer's machine (Ubuntu 24.04,
R 4.5.2, gcc 13), the two tarballs run back to back:

    checking tests    rejected tarball   73 s
    checking tests    this tarball       14 s

The same ratio holds for the suite run on its own: 90 s with `NOT_CRAN`
set against 16 s without it. Scaled onto the 12 minutes the rejected
tarball spent in `checking tests` on r-devel-windows-x86_64, this puts
the test block at roughly 2 minutes and the overall check comfortably
inside the 10-minute budget.

Coverage: with `NOT_CRAN=true` the suite runs 2522 tests, 3 skipped, all
passing. On CRAN it runs 2492 of them, 23 skipped.

## Test environments

* local: Ubuntu 24.04, R 4.5.2, gcc 13 (`R CMD check --as-cran`)

## R CMD check results

0 errors | 0 warnings | 0 notes attributable to the package.

The local run reports one warning and two notes, all three of which come
from the checking environment rather than from the package, and none of
which appeared on the CRAN pre-test machines:

* `'qpdf' is needed for checks on size reduction of PDFs` -- qpdf is not
  installed on the machine used for the check.
* `Skipping checking HTML validation: no command 'tidy' found` -- likewise.
* `Compilation used the following non-portable flag(s): '-mno-omit-leaf-frame-pointer'`
  -- this flag is injected by the Debian/Ubuntu build of R through
  `/usr/lib/R/etc/Makeconf`; it does not appear in the package.
  `src/Makevars` sets only `CXX_STD = CXX17`, `-I.`, and
  `$(SHLIB_OPENMP_CXXFLAGS)`.

## Compiled code

The package contains a C++17 kernel reached through Rcpp and
RcppArmadillo. OpenMP is requested through `$(SHLIB_OPENMP_CXXFLAGS)` so
that R resolves the correct flag for the toolchain, and every region
touching the runtime is guarded by `#ifdef _OPENMP`. The package builds
and runs single-threaded when OpenMP is unavailable, with identical
results.

## Package size

The previous development tarball carried a 6.6 MB tutorial video under
`man/figures/`. It has been excluded from the build; the tutorial is
linked from the README to its hosted location. The source tarball is
1.8 MB.

## Downstream dependencies

There are no reverse dependencies on CRAN.
