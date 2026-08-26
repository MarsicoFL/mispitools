## Submission of 2.0.1

This release answers the check failure reported for 2.0.0 on
r-release-macos-arm64, r-oldrel-macos-arm64 and the M1mac additional
check (23 test failures, all of them in the two files that cross-check
the C++ engine against the package's own R reference implementation).

### Diagnosis

The failures were not a difference in the computed likelihood ratios.
They were a difference in how many support points the distribution was
reported to have: the engine returned 37 203 atoms where the reference
returned 37 108, and the tests compare the two vectors elementwise.

An atom of the `log10` LR distribution is a real number that the two
implementations reach by different arithmetic routes. On aarch64 the
compiler contracts `a * b + c` into a fused multiply-add by default, so
one route rounds once where the other rounds twice, and two atoms that
are equal in exact arithmetic end up differing in their last bits.
Version 2.0.0 grouped atoms by IEEE equality, which then split one atom
into two. On x86-64, where the baseline ISA has no FMA and no
contraction takes place, both routes agree bit for bit and the tests
passed, which is why the failure was confined to the arm64 flavours.

The diagnosis was confirmed rather than inferred. Rebuilding 2.0.0 on
x86-64 with `-mfma -ffp-contract=fast` reproduces the reported failures
in the same two files and with the same signature (37 513 atoms against
37 251).

### Fix

Atom grouping now closes a group by a relative tolerance of 1e-12
instead of by exact equality, in the C++ kernel and in the R reference
alike, so the two agree on the size of the support on any platform. The
constant is chosen from a measurement: over the composed profiles the
suite exercises, consecutive keys are either within one unit in the last
place of each other, which means the same atom reached twice, or more
than 1e-8 apart in relative terms, which means genuinely distinct atoms.
The band between those two populations is empty across seven orders of
magnitude.

No compiler flag was added. `src/Makevars` still sets only
`CXX_STD = CXX17`, `-I.` and `$(SHLIB_OPENMP_CXXFLAGS)`.

The distribution itself is unchanged. Total mass is unchanged and the
expected weight of evidence agrees with 2.0.0 to twelve decimal places;
what changes is that duplicates split by rounding are now merged, so the
support is smaller and exact composition is correspondingly cheaper.

### Verification

Built and checked in two configurations on Ubuntu 24.04, R 4.5.2,
gcc 13: the ordinary build, and a build with `-mfma -mavx2
-ffp-contract=fast`, which stands in for the arm64 arithmetic that
produced the failures. Both pass the full suite with no failures and no
errors, in CRAN mode (2469 passing, 23 skipped) and with `NOT_CRAN=true`
(2519 passing, 3 skipped). Under 2.0.0 the same instrumented build
reproduces the CRAN failures.

Check time is unaffected by this release and remains well inside the
budget that 2.0.0 was resubmitted to meet; the exact composition path is
faster, since it now allocates a smaller support.

### On the short interval since 2.0.0

2.0.0 was published on 2026-08-25. This submission follows one day later
because it answers the check failure reported for it on the arm64
flavours, within the correction window given in the message of
2026-08-26. It contains that fix, its documentation, and nothing else.

## Test environments

* local: Ubuntu 24.04, R 4.5.2, gcc 13 (`R CMD check --as-cran`)
* local: the same, rebuilt with FMA contraction enabled to emulate arm64
* win-builder, R-devel

## R CMD check results

0 errors | 0 warnings | 0 notes attributable to the package.

The local run reports one warning and two notes, all three from the
checking environment rather than the package:

* `'qpdf' is needed for checks on size reduction of PDFs` -- qpdf is not
  installed on the machine used for the check.
* `Skipping checking HTML validation: no command 'tidy' found` -- likewise.
* `Compilation used the following non-portable flag(s): '-mno-omit-leaf-frame-pointer'`
  -- injected by the Debian/Ubuntu build of R through
  `/usr/lib/R/etc/Makeconf`; it does not appear in the package.

## Compiled code

A C++17 kernel reached through Rcpp and RcppArmadillo. OpenMP is
requested through `$(SHLIB_OPENMP_CXXFLAGS)` so that R resolves the
correct flag for the toolchain, and every region touching the runtime is
guarded by `#ifdef _OPENMP`. The package builds and runs single-threaded
when OpenMP is unavailable, with identical results.

## Downstream dependencies

There are no reverse dependencies on CRAN.
