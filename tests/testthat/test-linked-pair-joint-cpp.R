# F5.2 — core::linked_pair_joint() cross-checked against the canonical
# CRAN oracle pedprobr::likelihood2() (SCOUT_pedprobr_linkage.md). The
# roadmap names fbnet; SCOUT settled on pedprobr (cleaner, the official
# pedsuite dependency). Heavy mutation+linkage integral cross-checks are
# F5.3 (verifier); F5.2 validates the engine, the rho==0.5 factorisation,
# the map functions, the phase handling, and the boundary guards.

# --- test-local helpers ----------------------------------------------------

# Column-major-lex 1-based genotype index for an unordered allele pair,
# matching core/marker.h pair_to_genotype_index(): g(a1,a2) =
# a2*(a2+1)/2 + a1 with 0-based a1 <= a2, then +1.
.gidx <- function(i, j) {  # i, j are 1-based allele indices
  a1 <- min(i, j) - 1L
  a2 <- max(i, j) - 1L
  as.integer(a2 * (a2 + 1L) / 2L + a1 + 1L)
}

# Likelihood of an observed typed configuration under H1 (pedigree
# topology, recombination `rho`), obtained by marginalising the full
# linked-pair joint over the untyped members.
.linked_lik <- function(father, mother, freqs_a, freqs_b, rho, obs,
                        mut_a = list(0L, 0, 0, numeric(0)),
                        mut_b = list(0L, 0, 0, numeric(0))) {
  res <- cpp_linked_pair_joint(
    father = as.integer(father),
    mother = as.integer(mother),
    poi    = 0L,                       # H1 does not depend on the POI
    freqs_a = as.numeric(freqs_a),
    freqs_b = as.numeric(freqs_b),
    rho = rho,
    mutation_kind_a = mut_a[[1]], mutation_rate_a = mut_a[[2]],
    mutation_range_a = mut_a[[3]], numeric_labels_a = mut_a[[4]],
    mutation_kind_b = mut_b[[1]], mutation_rate_b = mut_b[[2]],
    mutation_range_b = mut_b[[3]], numeric_labels_b = mut_b[[4]])
  sa <- res$states_a
  sb <- res$states_b
  keep <- rep(TRUE, length(res$P_H1))
  for (nm in names(obs)) {
    mi <- as.integer(nm)              # 1-based member column
    g <- obs[[nm]]
    keep <- keep & (sa[, mi] == g[1]) & (sb[, mi] == g[2])
  }
  sum(res$P_H1[keep])
}

# Two-marker nuclear pedigree with `nch` children, all members typed,
# returns the pedprobr oracle likelihood at recombination `rho`.
.oracle_lik2 <- function(geno_a, geno_b, fa, fb, rho) {
  x <- pedtools::nuclearPed(length(geno_a) - 2L)
  ids <- as.character(seq_along(geno_a))
  x <- pedtools::addMarker(x, afreq = fa, name = "MA")
  x <- pedtools::addMarker(x, afreq = fb, name = "MB")
  for (k in seq_along(geno_a)) {
    pedtools::genotype(x, marker = "MA", id = ids[k]) <- geno_a[[k]]
    pedtools::genotype(x, marker = "MB", id = ids[k]) <- geno_b[[k]]
  }
  pedprobr::likelihood2(x, marker1 = 1, marker2 = 2, rho = rho)
}

# nuclearPed(n) flattened to 0-based father/mother (founders 1,2; the
# rest are their children). Member i (1-based) -> index i-1.
.nuclear_fm <- function(nch) {
  n <- nch + 2L
  father <- c(-1L, -1L, rep(0L, nch))
  mother <- c(-1L, -1L, rep(1L, nch))
  list(father = father, mother = mother, n = n)
}

# ---------------------------------------------------------------------------

test_that("map functions match pedprobr (haldane / kosambi)", {
  skip_if_not_installed("pedprobr")
  for (cM in c(1, 5, 12.5, 30, 80)) {
    expect_equal(haldane_cm_to_rho_wrap(cM),
                 pedprobr::haldane(cM = cM), tolerance = 1e-12)
    expect_equal(kosambi_cm_to_rho_wrap(cM),
                 pedprobr::kosambi(cM = cM), tolerance = 1e-12)
  }
  # Round trips.
  for (rho in c(0.01, 0.1, 0.3, 0.49)) {
    expect_equal(haldane_cm_to_rho_wrap(haldane_rho_to_cm_wrap(rho)),
                 rho, tolerance = 1e-12)
    expect_equal(kosambi_cm_to_rho_wrap(kosambi_rho_to_cm_wrap(rho)),
                 rho, tolerance = 1e-12)
  }
})

test_that("linked_pair_joint H1 == pedprobr::likelihood2 (rho sweep)", {
  skip_if_not_installed("pedtools")
  skip_if_not_installed("pedprobr")

  fa <- c("1" = 0.6, "2" = 0.4)
  fb <- c("1" = 0.7, "2" = 0.3)
  # nuclearPed(2): members 1,2 (parents) + 3,4 (children), all typed.
  geno_a <- list("1/2", "1/1", "1/2", "1/2")
  geno_b <- list("1/2", "1/1", "1/1", "1/2")
  fm <- .nuclear_fm(2L)

  # Observed genotype indices (allele "1" = index 1, "2" = index 2).
  obs <- list("1" = c(.gidx(1, 2), .gidx(1, 2)),
              "2" = c(.gidx(1, 1), .gidx(1, 1)),
              "3" = c(.gidx(1, 2), .gidx(1, 1)),
              "4" = c(.gidx(1, 2), .gidx(1, 2)))

  for (rho in c(0, 0.05, 0.1, 0.25, 0.5)) {
    ours <- .linked_lik(fm$father, fm$mother, fa, fb, rho, obs)
    oracle <- .oracle_lik2(geno_a, geno_b, fa, fb, rho)
    expect_equal(ours, oracle, tolerance = 1e-8,
                 info = paste0("rho = ", rho))
  }
})

test_that("rho == 0.5 factorises into the product of single-marker joints", {
  skip_if_not_installed("pedtools")
  skip_if_not_installed("pedprobr")

  fa <- c("1" = 0.55, "2" = 0.30, "3" = 0.15)
  fb <- c("1" = 0.5, "2" = 0.5)
  geno_a <- list("1/2", "2/3", "2/2")
  geno_b <- list("1/2", "1/1", "1/2")
  fm <- .nuclear_fm(1L)                       # trio

  obs <- list("1" = c(.gidx(1, 2), .gidx(1, 2)),
              "2" = c(.gidx(2, 3), .gidx(1, 1)),
              "3" = c(.gidx(2, 2), .gidx(1, 2)))

  ours <- .linked_lik(fm$father, fm$mother, fa, fb, 0.5, obs)
  # At rho = 0.5 the markers are unlinked: likelihood2 collapses to the
  # product of the two single-marker likelihoods.
  x1 <- pedtools::nuclearPed(1)
  x1 <- pedtools::addMarker(x1, `1` = "1/2", `2` = "2/3", `3` = "2/2",
                            afreq = fa, name = "MA")
  x2 <- pedtools::nuclearPed(1)
  x2 <- pedtools::addMarker(x2, `1` = "1/2", `2` = "1/1", `3` = "1/2",
                            afreq = fb, name = "MB")
  prod_single <- pedprobr::likelihood(x1, 1) * pedprobr::likelihood(x2, 1)
  expect_equal(ours, prod_single, tolerance = 1e-10)

  oracle <- .oracle_lik2(geno_a, geno_b, fa, fb, 0.5)
  expect_equal(ours, oracle, tolerance = 1e-10)
})

test_that("doubly-heterozygous founder phase is handled (vs pedprobr)", {
  skip_if_not_installed("pedtools")
  skip_if_not_installed("pedprobr")

  fa <- c("1" = 0.5, "2" = 0.5)
  fb <- c("1" = 0.5, "2" = 0.5)
  # Founder 1 is doubly heterozygous (1/2 at both loci): phase ambiguous.
  geno_a <- list("1/2", "1/2", "1/2")
  geno_b <- list("1/2", "1/2", "1/2")
  fm <- .nuclear_fm(1L)

  obs <- list("1" = c(.gidx(1, 2), .gidx(1, 2)),
              "2" = c(.gidx(1, 2), .gidx(1, 2)),
              "3" = c(.gidx(1, 2), .gidx(1, 2)))

  for (rho in c(0, 0.1, 0.25, 0.5)) {
    ours <- .linked_lik(fm$father, fm$mother, fa, fb, rho, obs)
    oracle <- .oracle_lik2(geno_a, geno_b, fa, fb, rho)
    expect_equal(ours, oracle, tolerance = 1e-8,
                 info = paste0("double-het rho = ", rho))
  }
})

test_that("rho == 0.5 mixed mutation factorises vs single-marker engine", {
  skip_if_not_installed("pedtools")

  ped <- pedtools::nuclearPed(1)
  fa <- c("12" = 0.2, "13" = 0.3, "14" = 0.5)   # numeric labels (stepwise)
  fb <- c("1" = 0.7, "2" = 0.3)
  fm <- .nuclear_fm(1L)

  # Single-marker reference joint (F2 engine), marker A stepwise.
  mod_a <- marker_model(ped, "MA", fa,
                        mutation = list(model = "stepwise",
                                        rate = 0.01, ratio = 0.3))
  jA <- cpt_marker_joint_cpp_wrap(mod_a, poi = "3")

  res <- cpp_linked_pair_joint(
    father = as.integer(fm$father), mother = as.integer(fm$mother),
    poi = 2L,
    freqs_a = as.numeric(fa), freqs_b = as.numeric(fb), rho = 0.5,
    mutation_kind_a = 2L, mutation_rate_a = 0.01,
    mutation_range_a = 0.3,
    numeric_labels_a = as.numeric(names(fa)),
    mutation_kind_b = 0L)

  # Marginal H1 of member 3 (child) at marker A from the linked joint
  # must equal the single-marker H1 marginal (rho = 0.5 ⇒ independent).
  mar_link <- tapply(res$P_H1, res$states_a[, 3], sum)
  mar_link <- mar_link / sum(mar_link)
  ga <- attr(jA, "alleles")
  jA_mar <- tapply(jA$P_H1, jA[["3"]], sum)
  jA_mar <- jA_mar / sum(jA_mar)
  # Reorder single-marker marginal onto genotype-index order.
  lab <- character(length(jA_mar))
  gi <- 0L
  for (j in seq_along(ga)) for (i in seq_len(j)) {
    gi <- gi + 1L; lab[gi] <- paste0(ga[i], "/", ga[j])
  }
  jA_vec <- as.numeric(jA_mar[lab])
  expect_equal(as.numeric(mar_link), jA_vec, tolerance = 1e-10)

  # Whole-joint sanity: P_H1 and P_H2 are proper distributions.
  expect_equal(sum(res$P_H1), 1, tolerance = 1e-10)
  expect_equal(sum(res$P_H2), 1, tolerance = 1e-10)
})

test_that("boundary guards: rho range, frequencies, pedigree", {
  fm <- .nuclear_fm(1L)
  expect_error(
    cpp_linked_pair_joint(fm$father, fm$mother, 0L,
                          c(0.5, 0.5), c(0.5, 0.5), rho = 0.7),
    "rho")
  expect_error(
    cpp_linked_pair_joint(fm$father, fm$mother, 0L,
                          c(0.5, 0.5), c(0.5, 0.5), rho = -0.01),
    "rho")
  expect_error(
    cpp_linked_pair_joint(fm$father, fm$mother, 0L,
                          c(0.5, 0.4), c(0.5, 0.5), rho = 0.1),
    "sum to 1")
  expect_error(
    cpp_linked_pair_joint(c(-1L, -1L), c(-1L), 0L,
                          c(0.5, 0.5), c(0.5, 0.5), rho = 0.1),
    "same length")
})
