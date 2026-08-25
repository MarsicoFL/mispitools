## F1.3 — Cross-validate the R reference engine `cpt_marker_joint_R()`
## (mutation = "none") against `pedprobr::oneMarkerDistribution()` on
## 3 canonical pedigrees x 5 Argentina markers.
##
## Oracle: pedprobr (peeling). Tolerance 1e-10 per the F1.3 spec.
## Observed agreement is at floating-point noise level (< 1e-18 in all
## 15 cases at the time of writing; tolerance is loose so that minor
## upstream changes in pedprobr or BLAS do not destabilize the gate).
##
## Argentina frequencies are trimmed to the top 4 most frequent alleles
## per marker via `top_k_freqs()` (helper-pedigrees.R). Trimming is a
## tractability bound for the brute-force R engine; correctness of the
## joint factorization is what is being verified, not coverage of the
## full Argentina allele table.

skip_if_no_oracle <- function() {
  testthat::skip_if_not_installed("pedtools")
  testthat::skip_if_not_installed("pedprobr")
}

joint_max_abs_diff <- function(ped, marker_name, freqs) {
  ids <- as.character(labels(ped))

  mm <- marker_model(ped, marker_name, freqs)
  tab <- mispitools:::cpt_marker_joint_R(mm)

  m_obj <- pedtools::marker(ped, afreq = freqs, name = marker_name)
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

ARG_MARKERS <- c("D10S1248", "D22S1045", "THO1", "D16S539", "D5S818")

PEDS <- list(
  parent_child           = ped_parent_child,
  half_sibs              = ped_half_sibs,
  grandparent_grandchild = ped_grandparent_grandchild
)

## Cost is set by pedigree size, not by the marker: at top_k = 4 the
## brute-force enumeration is 10^3 states for the 3-member parent_child
## and 10^5 for the two 5-member pedigrees, so the ten cells on those
## two carry ~95% of this file's runtime. The always-on subset keeps
## parent_child across all five markers (the marker axis) plus the
## half-sib cell on one marker (a topology with a non-trivial peeling
## order). The grandparent_grandchild topology is covered on CRAN by
## test-vs-pedprobr-cpp.R against the same oracle, so the remaining
## nine cells here are gated (see helper-cran.R).
CORE_MARKER <- "D10S1248"
CORE_LARGE_PED <- "half_sibs"

is_exhaustive_cell <- function(ped_name, marker_name) {
  if (ped_name == "parent_child") return(FALSE)
  !(ped_name == CORE_LARGE_PED && marker_name == CORE_MARKER)
}

for (ped_name in names(PEDS)) {
  for (marker_name in ARG_MARKERS) {
    local({
      pn <- ped_name
      mk <- marker_name
      test_that(sprintf("P_H1 matches pedprobr: %s / %s", pn, mk), {
        skip_if_no_oracle()
        if (is_exhaustive_cell(pn, mk)) skip_if_exhaustive_disabled()
        data(Argentina, package = "mispitools", envir = environment())
        ped <- PEDS[[pn]]()
        freqs <- top_k_freqs(Argentina, mk, top_k = 4L)
        d <- joint_max_abs_diff(ped, mk, freqs)
        expect_lt(d, 1e-10)
      })
    })
  }
}
