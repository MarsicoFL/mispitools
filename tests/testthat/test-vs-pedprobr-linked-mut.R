## F5.3 — Integral cross-check of the linked-pair engine with mutation
## AND linkage active simultaneously, against the canonical CRAN oracle
## pedprobr::likelihood2().
##
## Scope vs prior phases:
##   * F5.2 (test-linked-pair-joint-cpp.R) validated core::linked_pair_joint
##     against pedprobr::likelihood2 for nuclearPed topologies with EITHER
##     no mutation (rho sweep) OR mutation on a single marker at the
##     rho == 0.5 factorisation point.
##   * F5.3 closes the gap: two different mutation models on the two
##     linked markers at once (marker A stepwise, marker B equal-rate),
##     swept across rho, on a multi-child nuclear pedigree with a doubly
##     heterozygous transmitting parent so the likelihood is genuinely
##     recombination-informative (the oracle moves with rho).
##
## Oracle: pedprobr::likelihood2() with pedmut mutation matrices attached
## to each marker. The pedmut parametrisation mirrors the core exactly,
## established in F2.6 (test-vs-pedprobr-cpp.R):
##   equal    -> mutationMatrix("equal",    rate)
##   stepwise -> mutationMatrix("stepwise", rate, rate2 = 0, range = ratio)
##
## Tolerance: 1e-8 per ROADMAP "Verificación cruzada" (linkage row) and
## the F5.2 spec. Observed agreement at authoring time sat at or below
## 5.3e-23 in every cell (machine precision; the linked peeling and the
## pedprobr peeling accumulate in different orders, hence not exactly 0).
##
## NOTE ON THE NAMED first-cousin SCOPE: the STATE.md row for F5.3 names
## a first-cousin pedigree. core::linked_pair_joint() materialises the
## DENSE phase-conscious joint over *all* members of the pedigree (its
## return value is an n_rows x n_members table, n_rows growing as
## (genotypes_A * genotypes_B) ^ n_members). For the 8-member first-cousin
## pedigree this is ~9^8 ~ 4e7 rows even at K = 2 per marker and the call
## aborts with std::bad_alloc. This is an engine scalability limitation
## (the single-marker cpt engine has the same dense-joint design but was
## only ever exercised on <=5-member pedigrees). It is tracked as a new
## microtask in STATE.md (cpp_engine: peeled / typed-set marginalisation
## for linked_pair_joint) and an ESCALATION note. The first-cousin block
## below is written and skip()-ed so it activates once the engine can
## marginalise instead of densifying.

# --- file-local helpers (standalone per verifier rule 1) -------------------

# Column-major-lex 1-based genotype index for an unordered allele pair,
# matching core/marker.h pair_to_genotype_index().
.f53_gidx <- function(i, j) {
  a1 <- min(i, j) - 1L
  a2 <- max(i, j) - 1L
  as.integer(a2 * (a2 + 1L) / 2L + a1 + 1L)
}

# "a/b" with allele labels -> c(gidx_A, gidx_B) via per-marker label maps.
.f53_obs <- function(gA, gB, mapA, mapB) {
  to_idx <- function(s, mp) {
    p <- strsplit(s, "/", fixed = TRUE)[[1]]
    .f53_gidx(mp[[p[1]]], mp[[p[2]]])
  }
  obs <- list()
  for (id in names(gA)) {
    obs[[id]] <- c(to_idx(gA[[id]], mapA), to_idx(gB[[id]], mapB))
  }
  obs
}

# Likelihood of the observed two-marker configuration under H1
# (pedigree topology + recombination `rho`), obtained by marginalising
# the linked-pair joint over the untyped members. `relevant` is set to
# exactly the typed members so the F5.5 peeled engine sums every other
# member out instead of densifying the full joint (when every member is
# typed -- e.g. the nuclearPed cases -- this is the whole pedigree and
# the dense F5.2 path is recovered unchanged).
.f53_linked_lik <- function(father, mother, freqs_a, freqs_b, rho, obs,
                            mut_a, mut_b) {
  res <- cpp_linked_pair_joint(
    father = as.integer(father), mother = as.integer(mother),
    poi = 0L,                              # H1 does not depend on the POI
    freqs_a = as.numeric(freqs_a), freqs_b = as.numeric(freqs_b),
    rho = rho,
    mutation_kind_a = mut_a$kind, mutation_rate_a = mut_a$rate,
    mutation_range_a = mut_a$range,
    numeric_labels_a = mut_a$labels,
    mutation_kind_b = mut_b$kind, mutation_rate_b = mut_b$rate,
    mutation_range_b = mut_b$range,
    numeric_labels_b = mut_b$labels,
    relevant = as.integer(names(obs)))
  keep <- rep(TRUE, length(res$P_H1))
  for (nm in names(obs)) {
    mi <- as.integer(nm)
    g <- obs[[nm]]
    keep <- keep & (res$states_a[, mi] == g[1]) & (res$states_b[, mi] == g[2])
  }
  sum(res$P_H1[keep])
}

# pedtools ped carrying the two markers with pedmut mutation models and
# the typed genotypes; used for the pedprobr oracle.
.f53_oracle_ped <- function(ped, gA, gB, fa, fb, mutA, mutB) {
  x <- pedtools::addMarker(ped, afreq = fa, name = "MA", mutmod = mutA)
  x <- pedtools::addMarker(x,   afreq = fb, name = "MB", mutmod = mutB)
  for (id in names(gA)) {
    pedtools::genotype(x, marker = "MA", id = id) <- gA[[id]]
    pedtools::genotype(x, marker = "MB", id = id) <- gB[[id]]
  }
  x
}

# ---------------------------------------------------------------------------

test_that("nuclearPed(2): stepwise+equal mutation x linkage == pedprobr::likelihood2", {
  skip_if_not_installed("pedtools")
  skip_if_not_installed("pedprobr")
  skip_if_not_installed("pedmut")

  # nuclearPed(2): members 1,2 = parents, 3,4 = children.
  father <- c(-1L, -1L, 0L, 0L)
  mother <- c(-1L, -1L, 1L, 1L)

  # Marker A: integer-labelled STR, stepwise mutation (range matters at
  # K = 3). Marker B: K = 2, equal-rate mutation.
  fa <- c("12" = 0.5, "13" = 0.3, "14" = 0.2)
  fb <- c("1"  = 0.7, "2"  = 0.3)
  rateA <- 0.004; ratioA <- 0.1; rateB <- 0.005

  # Parent 1 is doubly heterozygous (1/2 at A *and* B); two children with
  # distinct genotypes -> the joint is recombination-informative, so the
  # oracle is a non-trivial function of rho (verified: it sweeps from
  # ~2.0e-8 at rho=0 to ~3.6e-7 at rho=0.5).
  gA <- list("1" = "12/13", "2" = "12/12", "3" = "12/13", "4" = "13/14")
  gB <- list("1" = "1/2",   "2" = "1/1",   "3" = "1/1",   "4" = "1/2")

  mapA <- c("12" = 1L, "13" = 2L, "14" = 3L)
  mapB <- c("1" = 1L, "2" = 2L)
  obs  <- .f53_obs(gA, gB, mapA, mapB)

  mutA_core <- list(kind = 2L, rate = rateA, range = ratioA,
                    labels = as.numeric(names(fa)))
  mutB_core <- list(kind = 1L, rate = rateB, range = 0,
                    labels = numeric(0))

  mutA <- pedmut::mutationMatrix("stepwise", alleles = names(fa),
                                 rate = rateA, rate2 = 0, range = ratioA)
  mutB <- pedmut::mutationMatrix("equal", alleles = names(fb),
                                 rate = rateB)
  x <- .f53_oracle_ped(pedtools::nuclearPed(2), gA, gB, fa, fb, mutA, mutB)

  ## Every rho costs one pedprobr::likelihood2() call on the oracle
  ## side, which is what makes this the slowest test in the file. The
  ## endpoints rho = 0 (complete linkage) and rho = 0.5 (independence)
  ## are the two that pin the model down; the interior points check
  ## monotonicity in between and are gated (see helper-cran.R).
  rhos <- if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    c(0, 0.1, 0.25, 0.5)
  } else {
    c(0, 0.5)
  }

  prev <- NA_real_
  for (rho in rhos) {
    ours   <- .f53_linked_lik(father, mother, fa, fb, rho, obs,
                              mutA_core, mutB_core)
    oracle <- pedprobr::likelihood2(x, marker1 = 1, marker2 = 2, rho = rho)
    expect_equal(ours, oracle, tolerance = 1e-8,
                 info = paste0("rho = ", rho))
    # Sanity: the configuration is genuinely linkage-informative, so the
    # likelihood is strictly increasing in rho here (not rho-invariant).
    if (!is.na(prev)) expect_gt(oracle, prev)
    prev <- oracle
  }
})

test_that("rho == 0.5 factorises into product of single-marker likelihoods (mutation on both)", {
  ## Redundant on CRAN: the rho = 0.5 cell of the test above already
  ## compares the same engine against pedprobr::likelihood2() on the
  ## same pedigree. This one restates the identity in factorised form,
  ## which is worth two more oracle calls locally but not on CRAN
  ## (see helper-cran.R).
  skip_if_exhaustive_disabled()
  skip_if_not_installed("pedtools")
  skip_if_not_installed("pedprobr")
  skip_if_not_installed("pedmut")

  father <- c(-1L, -1L, 0L, 0L)
  mother <- c(-1L, -1L, 1L, 1L)
  fa <- c("12" = 0.5, "13" = 0.3, "14" = 0.2)
  fb <- c("1"  = 0.7, "2"  = 0.3)
  rateA <- 0.004; ratioA <- 0.1; rateB <- 0.005
  gA <- list("1" = "12/13", "2" = "12/12", "3" = "12/13", "4" = "13/14")
  gB <- list("1" = "1/2",   "2" = "1/1",   "3" = "1/1",   "4" = "1/2")
  obs <- .f53_obs(gA, gB, c("12" = 1L, "13" = 2L, "14" = 3L),
                  c("1" = 1L, "2" = 2L))

  mutA <- pedmut::mutationMatrix("stepwise", alleles = names(fa),
                                 rate = rateA, rate2 = 0, range = ratioA)
  mutB <- pedmut::mutationMatrix("equal", alleles = names(fb),
                                 rate = rateB)

  ours <- .f53_linked_lik(
    father, mother, fa, fb, 0.5, obs,
    list(kind = 2L, rate = rateA, range = ratioA,
         labels = as.numeric(names(fa))),
    list(kind = 1L, rate = rateB, range = 0, labels = numeric(0)))

  xa <- pedtools::addMarker(pedtools::nuclearPed(2), afreq = fa,
                            name = "MA", mutmod = mutA)
  for (id in names(gA)) pedtools::genotype(xa, "MA", id) <- gA[[id]]
  xb <- pedtools::addMarker(pedtools::nuclearPed(2), afreq = fb,
                            name = "MB", mutmod = mutB)
  for (id in names(gB)) pedtools::genotype(xb, "MB", id) <- gB[[id]]

  prod_single <- pedprobr::likelihood(xa, 1) * pedprobr::likelihood(xb, 1)
  # 1e-9: at rho = 0.5 the two loci are unlinked, so likelihood2 must
  # collapse to the product of the two single-marker likelihoods even
  # with mutation on both. Observed |diff| ~ 5e-23.
  expect_equal(ours, prod_single, tolerance = 1e-9)
})

test_that("linked joint with mutation is a proper distribution (sums to 1)", {
  skip_if_not_installed("pedtools")

  father <- c(-1L, -1L, 0L, 0L)
  mother <- c(-1L, -1L, 1L, 1L)
  ## The joint has (K_a(K_a+1)/2 * K_b(K_b+1)/2)^n rows: 9^4 at
  ## K_a = K_b = 2, 18^4 at K_a = 3. Both exercise the same
  ## normalisation path, so the larger marker A runs off CRAN only
  ## (see helper-cran.R).
  fa <- if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    c("12" = 0.5, "13" = 0.3, "14" = 0.2)
  } else {
    c("12" = 0.6, "13" = 0.4)
  }
  fb <- c("1"  = 0.7, "2"  = 0.3)

  res <- cpp_linked_pair_joint(
    father, mother, poi = 3L,
    freqs_a = as.numeric(fa), freqs_b = as.numeric(fb), rho = 0.1,
    mutation_kind_a = 2L, mutation_rate_a = 0.004,
    mutation_range_a = 0.1,
    numeric_labels_a = as.numeric(names(fa)),
    mutation_kind_b = 1L, mutation_rate_b = 0.005,
    mutation_range_b = 0)
  expect_equal(sum(res$P_H1), 1, tolerance = 1e-9)
  expect_equal(sum(res$P_H2), 1, tolerance = 1e-9)
})

test_that("first-cousin x 2 linked markers (named F5.3 scope)", {
  ## Topology breadth (8 members) on top of the nuclearPed(2) cell that
  ## always runs; three more pedprobr::likelihood2() calls. Gated
  ## (see helper-cran.R).
  skip_if_exhaustive_disabled()
  # F5.5 closed the engine-scalability gap: core::linked_pair_joint()
  # now marginalises non-relevant members via Elston-Stewart online
  # variable elimination instead of materialising the dense joint over
  # all members. The 8-member first-cousin pedigree -- which aborted
  # with std::bad_alloc during F5.3 -- is verified here against the
  # canonical oracle pedprobr::likelihood2(), reactivating F5.3.
  skip_if_not_installed("pedtools")
  skip_if_not_installed("pedprobr")
  skip_if_not_installed("pedmut")

  # 8-member first cousins: GF(1) GM(2) SpouseA(3) SpouseB(4) Sib1(5)
  # Sib2(6) Cousin1(7) Cousin2(8). 0-based, founders = -1.
  father <- c(-1L, -1L, -1L, -1L, 0L, 0L, 2L, 3L)
  mother <- c(-1L, -1L, -1L, -1L, 1L, 1L, 4L, 5L)
  fa <- c("12" = 0.6, "13" = 0.4)
  fb <- c("1"  = 0.6, "2"  = 0.4)
  gA <- list("1" = "12/13", "7" = "12/12", "8" = "13/13")
  gB <- list("1" = "1/2",   "7" = "1/1",   "8" = "1/2")
  obs <- .f53_obs(gA, gB, c("12" = 1L, "13" = 2L), c("1" = 1L, "2" = 2L))

  mutA <- pedmut::mutationMatrix("stepwise", alleles = names(fa),
                                 rate = 0.004, rate2 = 0, range = 0.1)
  mutB <- pedmut::mutationMatrix("equal", alleles = names(fb),
                                 rate = 0.005)
  x <- pedtools::ped(id = 1:8,
                     fid = c(0, 0, 0, 0, 1, 1, 3, 4),
                     mid = c(0, 0, 0, 0, 2, 2, 5, 6),
                     sex = c(1, 2, 1, 1, 2, 2, 1, 1))
  x <- .f53_oracle_ped(x, gA, gB, fa, fb, mutA, mutB)

  ## As above: endpoints always, interior point off CRAN only.
  rhos <- if (identical(Sys.getenv("NOT_CRAN"), "true")) {
    c(0, 0.1, 0.5)
  } else {
    c(0, 0.5)
  }

  for (rho in rhos) {
    ours <- .f53_linked_lik(
      father, mother, fa, fb, rho, obs,
      list(kind = 2L, rate = 0.004, range = 0.1,
           labels = as.numeric(names(fa))),
      list(kind = 1L, rate = 0.005, range = 0, labels = numeric(0)))
    oracle <- pedprobr::likelihood2(x, marker1 = 1, marker2 = 2, rho = rho)
    expect_equal(ours, oracle, tolerance = 1e-8,
                 info = paste0("first-cousin rho = ", rho))
  }
})
