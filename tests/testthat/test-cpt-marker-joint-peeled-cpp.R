## F5.6 — peeled / typed-set-marginalising path for cpt_marker_joint_cpp().
##
## The default (relevant = empty) keeps the dense F2 engine bit-for-bit
## (covered by the other cpt tests). Here we exercise the new Elston-
## Stewart online variable-elimination path:
##   (a) peeled keep-all == dense full joint  (algorithm self-consistency)
##   (b) peeled subset    == dense full joint marginalised in R over the
##       same kept members  (exact marginal)
##   (c) POI-only marginal == HWE genotype prior  (analytic oracle)
##   (d) scales to a pedigree the dense engine cannot densify
##   (e) out-of-range `relevant` is rejected

## Aggregate a cpt_marker_joint_cpp() result to the marginal over a set
## of 0-based member columns: collapse identical kept-member genotype
## tuples, summing P_H1 / P_H2. Returns a data.frame sorted by the tuple.
.marginal_over <- function(res, members0) {
  st <- res$states
  if (is.null(dim(st))) st <- matrix(st, nrow = length(res$P_H1))
  cols <- members0 + 1L
  key <- apply(st[, cols, drop = FALSE], 1L, paste, collapse = ",")
  h1 <- tapply(res$P_H1, key, sum)
  h2 <- tapply(res$P_H2, key, sum)
  ord <- order(names(h1))
  data.frame(key = names(h1)[ord],
             P_H1 = as.numeric(h1)[ord],
             P_H2 = as.numeric(h2)[ord],
             stringsAsFactors = FALSE)
}

.expect_marginal_equal <- function(a, b, tol = 1e-12, info = NULL) {
  expect_identical(a$key, b$key, info = info)
  expect_equal(a$P_H1, b$P_H1, tolerance = tol, info = info)
  expect_equal(a$P_H2, b$P_H2, tolerance = tol, info = info)
}

## Three-generation pedigree (0-based, topologically ordered):
##   0 GF, 1 GM, 2 spouseA, 3 spouseB           founders
##   4 = child(0,1), 5 = child(0,1)             sibs
##   6 = child(4,2), 7 = child(5,3)             first cousins
.cousin_fm <- function() {
  list(father = c(-1L, -1L, -1L, -1L, 0L, 0L, 4L, 5L),
       mother = c(-1L, -1L, -1L, -1L, 1L, 1L, 2L, 3L))
}

test_that("peeled keep-all reproduces the dense full joint (K=2,3)", {
  fm <- .cousin_fm()
  n <- length(fm$father)
  ## The dense reference is 3^8 rows at K=2 and 6^8 = 1.7M rows at K=3;
  ## the K=3 pass is ~4x the cost of the whole rest of this file. K=2
  ## already exercises the keep-all peeling path, so K=3 is gated
  ## (see helper-cran.R).
  freq_sets <- list(c(a = 0.6, b = 0.4))
  if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    freq_sets <- c(freq_sets, list(c(a = 0.3, b = 0.5, c = 0.2)))
  }
  for (freqs in freq_sets) {
    dense <- cpt_marker_joint_cpp(fm$father, fm$mother, poi = 6L,
                                  freqs = unname(freqs))
    peeled <- cpt_marker_joint_cpp(fm$father, fm$mother, poi = 6L,
                                   freqs = unname(freqs),
                                   relevant = seq_len(n))
    all0 <- 0:(n - 1L)
    .expect_marginal_equal(.marginal_over(dense, all0),
                           .marginal_over(peeled, all0),
                           info = paste0("K=", length(freqs)))
  }
})

test_that("peeled subset == dense marginalised over the kept members", {
  fm <- .cousin_fm()
  freqs <- c(a = 0.4, b = 0.35, c = 0.25)        # K=3
  dense <- cpt_marker_joint_cpp(fm$father, fm$mother, poi = 6L,
                                freqs = unname(freqs))
  ## Typed = the two cousins (6, 7) + a grandparent (0); POI (6) is
  ## always forced into the keep-set by the engine.
  ## Each subset is checked against the same K=3 dense joint (6^8 rows).
  ## The two-cousin subset is the representative case; the other two
  ## vary the keep-set and are gated (see helper-cran.R).
  subsets <- list(c(6L, 7L))
  if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    subsets <- list(c(6L, 7L), c(0L, 6L, 7L), c(6L))
  }
  for (subset in subsets) {
    keep0 <- sort(unique(c(subset, 6L)))         # engine forces POI in
    peeled <- cpt_marker_joint_cpp(fm$father, fm$mother, poi = 6L,
                                   freqs = unname(freqs),
                                   relevant = subset + 1L)
    .expect_marginal_equal(
      .marginal_over(dense, keep0),
      .marginal_over(peeled, keep0),
      info = paste0("subset={", paste(subset, collapse = ","), "}"))
    ## Marginals are proper distributions.
    expect_equal(sum(peeled$P_H1), 1, tolerance = 1e-12)
    expect_equal(sum(peeled$P_H2), 1, tolerance = 1e-12)
  }
})

test_that("POI-only marginal equals the HWE genotype prior", {
  fm <- .cousin_fm()
  p <- 0.7; q <- 0.3
  res <- cpt_marker_joint_cpp(fm$father, fm$mother, poi = 6L,
                              freqs = c(p, q), relevant = 7L)  # 1-based POI
  m <- .marginal_over(res, 6L)
  ## K=2 genotypes in column-major-lex order: (0,0), (0,1), (1,1).
  hwe <- c(p * p, 2 * p * q, q * q)
  expect_equal(sort(m$P_H1), sort(hwe), tolerance = 1e-12)
  ## POI is a no-mutation descendant of HWE founders under H1 and an
  ## independent HWE founder under H2: both marginals are HWE.
  expect_equal(m$P_H1, m$P_H2, tolerance = 1e-12)
})

test_that("peeled scales where the dense joint is intractable", {
  ## 14-member, 4-allele pedigree: dense densifies G=10 over 14 members
  ## (10^14 rows). Peeled with two typed cousins stays trivial.
  father <- c(-1L, -1L, -1L, -1L, -1L, -1L,
              0L, 0L, 2L, 4L, 6L, 7L, 8L, 9L)
  mother <- c(-1L, -1L, -1L, -1L, -1L, -1L,
              1L, 1L, 3L, 5L, 1L, 3L, 5L, 1L)
  freqs <- c(0.4, 0.3, 0.2, 0.1)
  t0 <- proc.time()[["elapsed"]]
  res <- cpt_marker_joint_cpp(father, mother, poi = 12L,
                              freqs = freqs,
                              relevant = c(11L, 13L, 14L))  # 1-based
  dt <- proc.time()[["elapsed"]] - t0
  expect_lt(dt, 5)
  expect_equal(sum(res$P_H1), 1, tolerance = 1e-10)
  expect_equal(sum(res$P_H2), 1, tolerance = 1e-10)
  expect_true(all(is.finite(res$P_H1)))
  expect_true(all(res$P_H1 >= 0))
})

test_that("out-of-range relevant index is rejected", {
  fm <- .cousin_fm()
  expect_error(
    cpt_marker_joint_cpp(fm$father, fm$mother, poi = 6L,
                         freqs = c(0.6, 0.4), relevant = 99L),
    "relevant member index out of range")
})
