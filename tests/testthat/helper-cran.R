## Check-time budget on CRAN.
##
## CRAN's incoming auto-check rejects a submission whose overall check
## time exceeds 10 minutes on any flavour. The 2.0.0 submission of
## 2026-08-24 was rejected on exactly that count (17 min on
## r-devel-windows-x86_64), with the test suite accounting for 12 of
## those minutes.
##
## The cost is concentrated in the cross-engine verification files
## (test-vs-*.R). Those compare the mispitools engines against
## independent oracles -- pedprobr, pedmut, Familias, forrel -- over a
## full grid of pedigree x mutation model x marker. The slow cells are
## dominated by the oracle's brute-force enumeration over 5- to
## 8-member pedigrees, not by mispitools itself.
##
## No cell is deleted. Each grid is split in two: a representative
## subset that always runs, and the exhaustive remainder gated behind
## `skip_if_exhaustive_disabled()`. The remainder runs whenever
## `NOT_CRAN` is set to "true", which covers `devtools::test()`,
## `devtools::check()`, and CI. A developer therefore still sees the
## complete grid; CRAN sees a suite that fits the budget.
skip_if_exhaustive_disabled <- function() {
  testthat::skip_on_cran()
}

## Monte Carlo sample size for the forrel cross-checks: full precision
## off CRAN, a tenth of it on CRAN. Those gates are stated as multiples
## of the sample standard error of the mean, which is recomputed from
## the draws, so they remain correctly calibrated at the smaller N --
## what shrinks is the power to detect a small bias, not the validity
## of the comparison.
mc_sample_size <- function(full = 20000L, reduced = 2000L) {
  if (identical(Sys.getenv("NOT_CRAN"), "true")) full else reduced
}
