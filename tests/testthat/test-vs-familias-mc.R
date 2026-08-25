## F4.5 — Cross-validate the per-marker LR and the LR distribution
## against two independent forensic oracles:
##
##   Block A. `Familias::FamiliasPosterior(...)$LRperMarker` — the
##     scalar per-marker LR for a *fixed observed profile*. mispitools
##     produces the full distribution of the LR over the POI genotype;
##     evaluated at one fixed full profile g, P_H1(g) / P_H2(g) from
##     `cpt_marker_joint_cpp_wrap()` must equal the Familias LR for
##     that profile (the field-standard MP construction: H1 = POI
##     related per pedigree, H2 = POI an unrelated lone founder).
##     Tolerance 1e-8 (ROADMAP "Half-sib, mutation=equal" cell).
##
##   Block B. `forrel::profileSim()` — Monte Carlo. Profiles drawn
##     under H1 (the related pedigree, with the marker's mutation
##     model attached) are scored with mispitools' exact per-marker
##     LR; the empirical mean of log10 LR must match the exact
##     E[log10 LR | H1] = per_marker_kl()$e_log10_lr_h1 within a
##     statistical band (6 standard errors, fixed seed). A second
##     two-marker block checks that the Monte Carlo mean of the total
##     log10 LR matches the exact composed mean from
##     `lr_distribution()` (F4.4), which itself must equal the sum of
##     the per-marker means deterministically (1e-9).
##
## This closes the F4 verification chain: F1.3/F2.6 anchored the joint
## table to pedprobr; here the *LR* and its *distribution* are anchored
## to Familias (closed form) and forrel (Monte Carlo), the two oracles
## named in the F4.5 milestone. Both are Suggests-only.

skip_if_no_familias <- function() {
  testthat::skip_if_not_installed("pedtools")
  testthat::skip_if_not_installed("Familias")
}

skip_if_no_forrel <- function() {
  testthat::skip_if_not_installed("pedtools")
  testthat::skip_if_not_installed("forrel")
  testthat::skip_if_not_installed("pedmut")
}

## Canonicalize a "a/b" genotype string so its allele order matches the
## allele index order used by cpt_marker_joint_cpp_wrap() (= order of
## names(freqs)). Needed because forrel/Familias emit alleles in input
## order, not engine order.
canon_geno <- function(g, alleles) {
  a <- strsplit(g, "/", fixed = TRUE)[[1]]
  a <- a[order(match(a, alleles))]
  paste0(a[1L], "/", a[2L])
}

## P_H1(profile) / P_H2(profile) from the C++ joint table, evaluated at
## one fixed full profile (named character vector, names = member ids).
mispi_profile_lr <- function(model, profile) {
  tab <- mispitools:::cpt_marker_joint_cpp_wrap(model)
  alleles <- attr(tab, "alleles")
  members <- names(profile)
  prof_c <- vapply(profile, canon_geno, character(1), alleles = alleles)
  key_tab <- do.call(paste, c(tab[, members, drop = FALSE], list(sep = "|")))
  key_prof <- paste(prof_c, collapse = "|")
  idx <- which(key_tab == key_prof)
  testthat::expect_length(idx, 1L)
  tab$P_H1[idx] / tab$P_H2[idx]
}

## Familias per-marker LR for one fixed profile. `rel`/`unr` are
## FamiliasPedigree objects (related H1 / unrelated-POI H2); ref = 1
## (the unrelated pedigree), so LRperMarker[, "related"] = lik(H1)/lik(H2).
familias_profile_lr <- function(rel, unr, freqs, profile, mutation) {
  al <- names(freqs)
  if (mutation$model == "none") {
    loc <- Familias::FamiliasLocus(frequencies = as.numeric(freqs),
                                   name = "M", allelenames = al)
  } else if (mutation$model == "equal") {
    loc <- Familias::FamiliasLocus(
      frequencies = as.numeric(freqs), name = "M", allelenames = al,
      femaleMutationModel = "Equal", maleMutationModel = "Equal",
      femaleMutationRate = mutation$rate, maleMutationRate = mutation$rate)
  } else {
    stop("oracle covers mut none/equal only")
  }
  dm <- do.call(rbind, lapply(profile, function(g)
    strsplit(g, "/", fixed = TRUE)[[1]]))
  rownames(dm) <- names(profile)
  post <- Familias::FamiliasPosterior(
    list(unrelated = unr, related = rel), loc, dm, ref = 1)
  post$LRperMarker["M", "related"]
}

## ---- Canonical cases (pedtools ped + matching FamiliasPedigree) ----

## Trio (missing-person core): 1=father, 2=mother, 3=child=POI.
## H1: 3 child of 1,2.  H2: 3 unrelated lone founder (1,2 still founders).
trio_ped <- function() pedtools::nuclearPed(1)
trio_rel <- function() Familias::FamiliasPedigree(
  id = c("1", "2", "3"), dadid = c(NA, NA, "1"),
  momid = c(NA, NA, "2"), sex = c("male", "female", "male"))
trio_unr <- function() Familias::FamiliasPedigree(
  id = c("1", "2", "3"), dadid = c(NA, NA, NA),
  momid = c(NA, NA, NA), sex = c("male", "female", "male"))

## Full sibs: 1,2 parents of 3 and 4; POI = 4 (engine default).
## H2: 3 stays a child of 1,2; 4 becomes an unrelated lone founder,
## mirroring P_H2 = P_H1(g1,g2,g3) * P_HWE(g4).
fs_ped <- function() pedtools::nuclearPed(2)
fs_rel <- function() Familias::FamiliasPedigree(
  id = c("1", "2", "3", "4"), dadid = c(NA, NA, "1", "1"),
  momid = c(NA, NA, "2", "2"), sex = c("male", "female", "male", "male"))
fs_unr <- function() Familias::FamiliasPedigree(
  id = c("1", "2", "3", "4"), dadid = c(NA, NA, "1", NA),
  momid = c(NA, NA, "2", NA), sex = c("male", "female", "male", "male"))

## Hardcoded frequencies (reproducibility; verifier rule 2). Integer
## allele labels so the stepwise path could reuse them if extended.
FR3 <- c("12" = 0.2, "13" = 0.3, "14" = 0.5)
FR4 <- c("15" = 0.1, "16" = 0.25, "17" = 0.3, "18" = 0.35)

## Mendelian-consistent profiles (so the LR is finite under mut="none":
## P_H1 > 0). For the trio POI=3, for full-sibs POI=4.
TRIO_PROF3 <- c("1" = "12/13", "2" = "13/14", "3" = "12/14")
TRIO_PROF4 <- c("1" = "15/16", "2" = "16/17", "3" = "15/17")
FS_PROF3   <- c("1" = "12/13", "2" = "13/14", "3" = "12/13", "4" = "13/14")
FS_PROF4   <- c("1" = "15/16", "2" = "16/17", "3" = "15/16", "4" = "16/17")

MUT_NONE  <- list(model = "none", rate = 0)
MUT_EQUAL <- list(model = "equal", rate = 0.004)

FAM_CASES <- list(
  list(nm = "trio / FR3", pedf = trio_ped, relf = trio_rel,
       unrf = trio_unr, fr = FR3, prof = TRIO_PROF3),
  list(nm = "trio / FR4", pedf = trio_ped, relf = trio_rel,
       unrf = trio_unr, fr = FR4, prof = TRIO_PROF4),
  list(nm = "fullsib / FR3", pedf = fs_ped, relf = fs_rel,
       unrf = fs_unr, fr = FR3, prof = FS_PROF3),
  list(nm = "fullsib / FR4", pedf = fs_ped, relf = fs_rel,
       unrf = fs_unr, fr = FR4, prof = FS_PROF4)
)

for (cs in FAM_CASES) {
  for (mut_nm in c("none", "equal")) {
    local({
      cs_ <- cs
      mut_ <- if (mut_nm == "none") MUT_NONE else MUT_EQUAL
      label <- sprintf("%s / mut=%s", cs_$nm, mut_$model)
      test_that(sprintf("Familias LRperMarker == P_H1/P_H2: %s", label), {
        skip_if_no_familias()
        ped <- cs_$pedf()
        mm <- marker_model(ped, "M", cs_$fr, mutation = mut_)
        lr_mispi <- mispi_profile_lr(mm, cs_$prof)
        lr_fam <- familias_profile_lr(cs_$relf(), cs_$unrf(),
                                      cs_$fr, cs_$prof, mut_)
        ## Tolerance 1e-8 per ROADMAP. Observed agreement is at
        ## floating-point noise (<= 1e-12) because both reduce to the
        ## same Elston-Stewart factorization on a 3-4 member pedigree.
        expect_equal(lr_mispi, unname(lr_fam), tolerance = 1e-8)
      })
    })
  }
}

## ---- Block B: Monte Carlo (forrel::profileSim) ----

## Empirical mean of log10 LR over N profiles simulated under H1.
mc_mean_log10lr <- function(ped, freqs, mutation, n, seed) {
  members <- as.character(labels(ped))
  mm <- marker_model(ped, "M", freqs, mutation = mutation)
  tab <- mispitools:::cpt_marker_joint_cpp_wrap(mm)
  alleles <- attr(tab, "alleles")
  key_tab <- do.call(paste, c(tab[, members, drop = FALSE],
                              list(sep = "|")))
  lr_lookup <- stats::setNames(tab$P_H1 / tab$P_H2, key_tab)

  mk <- if (mutation$model == "equal") {
    mut <- pedmut::mutationMatrix("equal", alleles = names(freqs),
                                  rate = mutation$rate)
    pedtools::marker(ped, afreq = freqs, name = "M", mutmod = mut)
  } else {
    pedtools::marker(ped, afreq = freqs, name = "M")
  }
  pedM <- pedtools::addMarkers(ped, mk)

  set.seed(seed)
  sims <- suppressMessages(forrel::profileSim(
    pedM, N = n, ids = members, verbose = FALSE))
  ll <- vapply(sims, function(s) {
    gg <- pedtools::getGenotypes(s, ids = members)[, 1]
    gg <- vapply(gg, canon_geno, character(1), alleles = alleles)
    log10(lr_lookup[[paste(gg, collapse = "|")]])
  }, numeric(1))
  list(mean = mean(ll), se = stats::sd(ll) / sqrt(n))
}

test_that("forrel Monte Carlo: E[log10 LR|H1] matches per_marker_kl (mut=none)", {
  skip_if_no_forrel()
  ped <- pedtools::nuclearPed(1)
  mm <- marker_model(ped, "M", FR3, mutation = MUT_NONE)
  exact <- per_marker_kl(mm)$e_log10_lr_h1
  mc <- mc_mean_log10lr(ped, FR3, MUT_NONE, n = mc_sample_size(), seed = 4071L)
  ## Statistical tolerance: 6 standard errors of the sample mean. The
  ## sampler (profileSim) is independent of the exact engine; this is a
  ## distribution-shape check, not a deterministic identity. 6*SE keeps
  ## the gate from flaking while still failing on real bias.
  expect_lt(abs(mc$mean - exact), 6 * mc$se)
})

test_that("forrel Monte Carlo: E[log10 LR|H1] matches per_marker_kl (mut=equal)", {
  skip_if_no_forrel()
  ped <- pedtools::nuclearPed(1)
  mm <- marker_model(ped, "M", FR3, mutation = MUT_EQUAL)
  exact <- per_marker_kl(mm)$e_log10_lr_h1
  mc <- mc_mean_log10lr(ped, FR3, MUT_EQUAL, n = mc_sample_size(), seed = 4072L)
  expect_lt(abs(mc$mean - exact), 6 * mc$se)
})

test_that("forrel Monte Carlo: composed two-marker mean matches lr_distribution()", {
  skip_if_no_forrel()
  ped <- pedtools::nuclearPed(1)
  members <- as.character(labels(ped))

  m1 <- marker_model(ped, "MA", FR3, mutation = MUT_EQUAL)
  m2 <- marker_model(ped, "MB", FR4, mutation = MUT_NONE)

  ## Deterministic anchor: the composed exact mean (F4.4 lr_distribution
  ## -> F4.2 convolution) equals the sum of the per-marker exact means
  ## (E is additive under conditional independence). 1e-9.
  exact_total <- per_marker_kl(m1)$e_log10_lr_h1 +
    per_marker_kl(m2)$e_log10_lr_h1
  comp_mean <- summary(lr_distribution(list(m1, m2)))$mean_h1
  expect_equal(comp_mean, exact_total, tolerance = 1e-9)

  ## Monte Carlo: simulate both markers jointly under H1, sum the
  ## per-marker log10 LR per profile, compare empirical mean to the
  ## exact composed mean within 6 SE.
  tab1 <- mispitools:::cpt_marker_joint_cpp_wrap(m1)
  tab2 <- mispitools:::cpt_marker_joint_cpp_wrap(m2)
  al1 <- attr(tab1, "alleles"); al2 <- attr(tab2, "alleles")
  lk1 <- stats::setNames(
    tab1$P_H1 / tab1$P_H2,
    do.call(paste, c(tab1[, members, drop = FALSE], list(sep = "|"))))
  lk2 <- stats::setNames(
    tab2$P_H1 / tab2$P_H2,
    do.call(paste, c(tab2[, members, drop = FALSE], list(sep = "|"))))

  mutA <- pedmut::mutationMatrix("equal", alleles = names(FR3),
                                 rate = MUT_EQUAL$rate)
  pedM <- pedtools::addMarkers(
    ped, pedtools::marker(ped, afreq = FR3, name = "MA", mutmod = mutA))
  pedM <- pedtools::addMarkers(
    pedM, pedtools::marker(ped, afreq = FR4, name = "MB"))

  set.seed(4073L)
  n <- mc_sample_size()
  sims <- suppressMessages(forrel::profileSim(
    pedM, N = n, ids = members, verbose = FALSE))
  tot <- vapply(sims, function(s) {
    gg <- pedtools::getGenotypes(s, ids = members)
    g1 <- vapply(gg[, "MA"], canon_geno, character(1), alleles = al1)
    g2 <- vapply(gg[, "MB"], canon_geno, character(1), alleles = al2)
    log10(lk1[[paste(g1, collapse = "|")]]) +
      log10(lk2[[paste(g2, collapse = "|")]])
  }, numeric(1))
  se <- stats::sd(tot) / sqrt(n)
  expect_lt(abs(mean(tot) - exact_total), 6 * se)
})
