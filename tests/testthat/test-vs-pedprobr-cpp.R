## F2.6 — Exhaustive cross-check of cpt_marker_joint_cpp_wrap() against
## an independent pedprobr + pedmut oracle. 12 cells = 3 pedigrees x
## 2 mutation models x 2 marker sets.
##
## Scope vs prior phases:
##   * F1.3 validated cpt_marker_joint_R (mut="none") vs pedprobr.
##   * F1.5 validated cpt_marker_joint_R (mut="stepwise") vs pedprobr+pedmut.
##   * F2.2 / F2.4 validated cpt_marker_joint_cpp vs the R reference
##     engine bit-for-bit (1e-12) across mut=none/equal/stepwise.
## F2.6 closes the chain by hitting the C++ pipeline (Rcpp binding +
## src/core peeling + mutation matrix builders) directly with a third
## independent peeling implementation. This catches a binding-layer or
## mutation-builder bug that would otherwise hide behind a shared
## R-ref oracle.
##
## Tolerance: 1e-10 per the F1.3 spec. Observed agreement sits near
## floating-point noise (<1e-14 in all 12 cells at the time of writing).
##
## Argentina frequencies are trimmed by top_k_freqs() (helper-pedigrees.R)
## so that the brute-force enumeration that backs the comparator stays
## tractable. Trimming bounds state space; it does not weaken the
## correctness check, since the joint factorization is what is being
## validated, not allele coverage.

skip_if_no_oracle <- function() {
  testthat::skip_if_not_installed("pedtools")
  testthat::skip_if_not_installed("pedprobr")
  testthat::skip_if_not_installed("pedmut")
}

build_pedmut_matrix <- function(alleles, mut_spec) {
  switch(mut_spec$model,
    equal = pedmut::mutationMatrix(
      model = "equal", alleles = alleles, rate = mut_spec$rate
    ),
    stepwise = pedmut::mutationMatrix(
      model = "stepwise", alleles = alleles,
      rate = mut_spec$rate, rate2 = 0, range = mut_spec$ratio
    ),
    stop("F2.6 oracle covers mut=equal and mut=stepwise only.")
  )
}

joint_max_abs_diff_cpp <- function(ped, marker_name, freqs, mut_spec) {
  ids <- as.character(labels(ped))

  mm <- marker_model(ped, marker_name, freqs, mutation = mut_spec)
  tab <- mispitools:::cpt_marker_joint_cpp_wrap(mm)

  pedmut_mat <- build_pedmut_matrix(names(freqs), mut_spec)
  m_obj <- pedtools::marker(ped, afreq = freqs, name = marker_name,
                            mutmod = pedmut_mat)
  ped_with <- pedtools::addMarkers(ped, m_obj)
  oracle <- pedprobr::oneMarkerDistribution(
    ped_with, ids = ids, marker = marker_name,
    output = "table", verbose = FALSE
  )

  key_m <- do.call(paste, c(tab[, ids, drop = FALSE], list(sep = "|")))
  key_o <- do.call(paste, c(oracle[, ids, drop = FALSE], list(sep = "|")))
  m_prob <- tapply(tab$P_H1, key_m, sum)
  o_prob <- stats::setNames(oracle$prob, key_o)

  all_keys <- union(names(m_prob), names(o_prob))
  m_vec <- unname(m_prob[all_keys]); m_vec[is.na(m_vec)] <- 0
  o_vec <- unname(o_prob[all_keys]); o_vec[is.na(o_vec)] <- 0
  max(abs(m_vec - o_vec))
}

PEDS_F26 <- list(
  parent_child           = ped_parent_child,
  half_sibs              = ped_half_sibs,
  grandparent_grandchild = ped_grandparent_grandchild
)

MUTS_F26 <- list(
  equal    = list(model = "equal",    rate = 0.005),
  stepwise = list(model = "stepwise", rate = 0.005, ratio = 0.1)
)

## Two marker sets pulled from Argentina with different K to exercise
## the joint table size axis. D10S1248 -> top-3 (K=3, 6 genotypes per
## member). D5S818 -> top-4 (K=4, 10 genotypes per member). Both are
## STR loci with integer-coded alleles, valid for the stepwise model.
##
## The top-4 set is the expensive half of the grid: the pedprobr oracle
## enumerates 10 genotypes per member over 5-member pedigrees, and those
## six cells alone account for ~95% of this file's runtime. They are
## therefore marked `exhaustive` and gated (see helper-cran.R). The
## top-3 set keeps all three pedigrees and both mutation models covered
## on CRAN, which is what closes the binding-layer chain described above.
MARKER_SETS_F26 <- list(
  "D10S1248_top3" = list(marker = "D10S1248", top_k = 3L,
                         exhaustive = FALSE),
  "D5S818_top4"   = list(marker = "D5S818",   top_k = 4L,
                         exhaustive = TRUE)
)

for (pn in names(PEDS_F26)) {
  for (mut_name in names(MUTS_F26)) {
    for (ms_name in names(MARKER_SETS_F26)) {
      local({
        pn_ <- pn; mut_ <- mut_name; ms_ <- ms_name
        test_that(sprintf("cpp vs pedprobr+pedmut: %s / %s / %s",
                          pn_, mut_, ms_), {
          skip_if_no_oracle()
          ms <- MARKER_SETS_F26[[ms_]]
          if (ms$exhaustive) skip_if_exhaustive_disabled()
          data(Argentina, package = "mispitools", envir = environment())
          ped <- PEDS_F26[[pn_]]()
          freqs <- top_k_freqs(Argentina, ms$marker, top_k = ms$top_k)
          d <- joint_max_abs_diff_cpp(ped, ms$marker, freqs,
                                      MUTS_F26[[mut_]])
          expect_lt(d, 1e-10)
        })
      })
    }
  }
}
