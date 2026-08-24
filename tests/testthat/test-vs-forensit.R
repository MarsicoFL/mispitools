## F1.7 — Cross-validate `per_marker_kl_R()` against the "forensIT-style"
## per-marker MP-only marginal KL of the missing person's genotype.
##
## Oracle (verifier note, 2026-05-13):
##
## The F1.7 spec named `forensIT::perMarkerKLs()` as the oracle.
## forensIT v1.2.0 has a defect in that wrapper: in the loop body
##
##     pop <- as.data.frame(frequency[i])
##     pop$Allele <- rownames(pop)
##     names(pop) <- c("Allele", "freq")
##
## the `names(pop)` reassignment is positional, so after the rename
## column 1 ("Allele") holds the numeric frequency values and column 2
## ("freq") holds the allele labels (originally rownames). The
## subsequent `pop$freq[match(Allele1, pop$Allele)]` lookup never
## resolves; every RPT entry becomes NA; and `replace(df, is.na(df) |
## df == 0, 1e-20)` collapses RPT to 1e-20 everywhere. The function
## then returns ~20 for any deterministic CPT, regardless of the true
## KL. The bundled example in `?perMarkerKLs` is itself affected
## (KLpedpop ≈ 18-19 for the Norwegian-frequency linearPed(2)
## profile).
##
## The active comparison below targets the *math* perMarkerKLs was
## meant to implement, by rebuilding the wrapper from the same
## underlying pedprobr machinery it uses internally:
##
##     CPT(g_MP) = pedprobr::oneMarkerDistribution(typed_ped, MP)
##     RPT(g_MP) = HWE under the population allele frequencies
##     KL_log10  = sum_{g_MP : CPT > 0} CPT * (log10(CPT) - log10(RPT))
##
## The exact chain-rule identity verified here (cf. F1.6 e_log10_lr_h1
## construction) is:
##
##     joint_KL(P_H1 || P_H2) in log10 base
##         = E_{typings of all non-MP members ~ P_H1}
##           [ KL( P_H1(g_MP | typings) || HWE(g_MP) ) in log10 ]
##
## valid whenever MP is a leaf node (so P_H1(typings) = P_H2(typings)
## and the typing-margin term drops out of the joint KL). All MP IDs
## chosen below are leaves of the canonical pedigrees.
##
## Tolerance: 1e-10 (the F1.7 spec is 1e-8; observed agreement is
## floating-point ULP-level, < 5e-16 in every case checked at the time
## of writing). The looser spec tolerance accommodates upstream
## pedprobr / BLAS perturbations.
##
## A separate `test_that` block exercises forensIT::perMarkerKLs()
## directly; it is `skip()`-ed pending the upstream fix. Removing the
## skip in a future forensIT will yield a passing test.

skip_if_no_deps <- function() {
  testthat::skip_if_not_installed("pedtools")
  testthat::skip_if_not_installed("pedprobr")
  testthat::skip_if_not_installed("pedmut")
}

skip_if_no_forensit <- function() {
  testthat::skip_if_not_installed("forensIT")
}

## -- helpers ----------------------------------------------------------

# All ordered genotype labels and their HWE weights for `freqs`.
build_geno_table <- function(freqs) {
  a <- names(freqs)
  geno_list <- list()
  hwe_vec <- numeric(0)
  geno_labels <- character(0)
  for (i in seq_along(a)) {
    for (j in i:length(a)) {
      geno_list[[length(geno_list) + 1L]] <- c(a[i], a[j])
      hwe_vec <- c(hwe_vec, if (i == j) freqs[i]^2 else 2 * freqs[i] * freqs[j])
      geno_labels <- c(geno_labels, paste0(a[i], "/", a[j]))
    }
  }
  names(hwe_vec) <- geno_labels
  list(geno_list = geno_list, hwe_vec = hwe_vec, geno_labels = geno_labels)
}

# Build a marker with optional mutation model and set it on `ped`.
typed_ped <- function(ped, typings, freqs, mutmod = NULL) {
  args <- c(list(ped, afreq = freqs, name = "TM"), typings)
  if (!is.null(mutmod)) args$mutmod <- mutmod
  m <- do.call(pedtools::marker, args)
  pedtools::setMarkers(ped, m)
}

# Chain-rule "forensIT-equivalent" weighted KL:
#   joint_KL_log10 = E_{typings of all non-MP members ~ P_H1}
#                      [ KL( P_H1(g_MP | typings) || HWE(g_MP) ) in log10 ]
# Weights are computed via pedprobr::likelihood() on the typed ped.
# When `directions = "both"`, also returns the H2->H1 average (only
# finite when CPT > 0 everywhere RPT > 0; use with mutation > 0).
forensit_equivalent_kl <- function(ped, mp_id, freqs, mutmod = NULL,
                                    directions = c("h1", "both")) {
  directions <- match.arg(directions)
  members <- as.character(labels(ped))
  non_mp <- setdiff(members, mp_id)
  gt <- build_geno_table(freqs)
  G <- length(gt$geno_list)
  grid <- as.matrix(do.call(expand.grid,
                             replicate(length(non_mp), seq_len(G),
                                       simplify = FALSE)))
  colnames(grid) <- non_mp
  acc_h1 <- 0; acc_h2 <- 0
  for (r in seq_len(nrow(grid))) {
    typings <- lapply(non_mp, function(id) gt$geno_list[[grid[r, id]]])
    names(typings) <- non_mp
    x <- typed_ped(ped, typings, freqs, mutmod = mutmod)
    w <- pedprobr::likelihood(x, marker = 1, verbose = FALSE)
    if (!is.finite(w) || w <= 0) next
    cpt <- pedprobr::oneMarkerDistribution(x, ids = mp_id, marker = 1,
                                            verbose = FALSE)
    cpt_v <- as.numeric(cpt); names(cpt_v) <- rownames(cpt)
    rpt_v <- gt$hwe_vec[names(cpt_v)]
    kl_h1 <- sum(ifelse(cpt_v > 0,
                         cpt_v * (log10(cpt_v) - log10(rpt_v)), 0))
    acc_h1 <- acc_h1 + w * kl_h1
    if (directions == "both") {
      kl_h2 <- sum(ifelse(rpt_v > 0 & cpt_v > 0,
                           rpt_v * (log10(rpt_v) - log10(cpt_v)), 0))
      acc_h2 <- acc_h2 + w * kl_h2
    }
  }
  if (directions == "h1") {
    c(kl_h1_log10 = acc_h1)
  } else {
    c(kl_h1_log10 = acc_h1, kl_h2_log10 = acc_h2)
  }
}

## -- tests ------------------------------------------------------------

test_that("parent_child trio, 2 alleles, no mutation: e_log10_lr_h1 matches forensIT-equivalent oracle", {
  # Tolerance: 1e-10. Observed: ~5e-17 (1 ULP). Forensic interpretation:
  # average KL of the child's marginal given parents typed exhaustively
  # under HWE, against the HWE marginal — the canonical MP-trio
  # information measure.
  skip_if_no_deps()
  freqs <- c("1" = 0.4, "2" = 0.6)
  ped <- pedtools::nuclearPed(1)
  mp <- "3"
  mm <- mispitools::marker_model(ped, "TM", freqs)
  my_kl <- mispitools:::per_marker_kl_R(mm, poi = mp)
  oracle <- forensit_equivalent_kl(ped, mp_id = mp, freqs = freqs)
  expect_equal(my_kl$e_log10_lr_h1, unname(oracle["kl_h1_log10"]),
               tolerance = 1e-10)
})

test_that("parent_child trio, 3 alleles, no mutation: e_log10_lr_h1 matches forensIT-equivalent oracle", {
  # 3-allele variant of the trio. Enumerated typings = 6^2 = 36.
  skip_if_no_deps()
  freqs <- c("1" = 0.3, "2" = 0.5, "3" = 0.2)
  ped <- pedtools::nuclearPed(1)
  mp <- "3"
  mm <- mispitools::marker_model(ped, "TM", freqs)
  my_kl <- mispitools:::per_marker_kl_R(mm, poi = mp)
  oracle <- forensit_equivalent_kl(ped, mp_id = mp, freqs = freqs)
  expect_equal(my_kl$e_log10_lr_h1, unname(oracle["kl_h1_log10"]),
               tolerance = 1e-10)
})

test_that("half_sibs, 2 alleles, no mutation: e_log10_lr_h1 matches forensIT-equivalent oracle", {
  # halfSibPed: members 1..5; founders 1,2,3; nonfounders 4 (child of
  # 1,2) and 5 (child of 2,3). MP = 5 (leaf). Enumerated typings on
  # the other 4 members = 3^4 = 81.
  skip_if_no_deps()
  freqs <- c("1" = 0.4, "2" = 0.6)
  ped <- pedtools::halfSibPed()
  mp <- "5"
  mm <- mispitools::marker_model(ped, "TM", freqs)
  my_kl <- mispitools:::per_marker_kl_R(mm, poi = mp)
  oracle <- forensit_equivalent_kl(ped, mp_id = mp, freqs = freqs)
  expect_equal(my_kl$e_log10_lr_h1, unname(oracle["kl_h1_log10"]),
               tolerance = 1e-10)
})

test_that("grandparent_grandchild (linearPed(2)), 2 alleles, no mutation: e_log10_lr_h1 matches forensIT-equivalent oracle", {
  # linearPed(2): members 1..5; founders 1,2,4; nonfounders 3 (child
  # of 1,2) and 5 (child of 3,4). MP = 5 (leaf). The intermediate
  # nonfounder 3 must be typed and its Mendelian-conditional weight
  # included in the chain-rule average — this is the case where the
  # "type only founders" shortcut fails (data-processing inequality).
  # 3^4 = 81 typings with many Mendelian-zero weight.
  skip_if_no_deps()
  freqs <- c("1" = 0.4, "2" = 0.6)
  ped <- pedtools::linearPed(2)
  mp <- "5"
  mm <- mispitools::marker_model(ped, "TM", freqs)
  my_kl <- mispitools:::per_marker_kl_R(mm, poi = mp)
  oracle <- forensit_equivalent_kl(ped, mp_id = mp, freqs = freqs)
  expect_equal(my_kl$e_log10_lr_h1, unname(oracle["kl_h1_log10"]),
               tolerance = 1e-10)
})

test_that("parent_child trio, 3 alleles, mutation=equal rate=0.005: both KL directions match", {
  # With strictly positive mutation, P_H1(g_c | g_f, g_m) > 0 for every
  # g_c, so the H2->H1 direction is finite and verifiable too.
  # Oracle uses `pedmut::mutationModel("equal", ...)`, the same
  # construction as the F1.4 equal-rate mutation_matrix_R().
  skip_if_no_deps()
  freqs <- c("1" = 0.3, "2" = 0.5, "3" = 0.2)
  ped <- pedtools::nuclearPed(1)
  mp <- "3"
  rate <- 0.005
  mut_mat <- pedmut::mutationModel("equal", alleles = names(freqs),
                                     rate = rate)
  mm <- mispitools::marker_model(ped, "TM", freqs,
                                  mutation = list(model = "equal",
                                                  rate = rate))
  my_kl <- mispitools:::per_marker_kl_R(mm, poi = mp)
  oracle <- forensit_equivalent_kl(ped, mp_id = mp, freqs = freqs,
                                    mutmod = mut_mat, directions = "both")
  expect_equal(my_kl$e_log10_lr_h1, unname(oracle["kl_h1_log10"]),
               tolerance = 1e-10)
  # e_log10_lr_h2 = E_H2[log10 LR] = -KL(H2 || H1) in log10 base, so
  # oracle's "kl_h2_log10" (= KL(RPT||CPT) in log10) equals
  # -my$e_log10_lr_h2.
  expect_equal(my_kl$e_log10_lr_h2, -unname(oracle["kl_h2_log10"]),
               tolerance = 1e-10)
})

test_that("documentation probe: forensIT::perMarkerKLs (skipped pending v1.2.0 RPT lookup fix)", {
  # Direct invocation of the F1.7-spec oracle. Skipped because
  # forensIT v1.2.0 mis-renames the frequency data.frame columns,
  # collapsing every RPT entry to 1e-20 and returning ~20 for any
  # deterministic CPT. See header for full diagnosis. When upstream
  # forensIT is patched, drop the skip() below; the comparison should
  # then agree with the active tests at 1e-10.
  skip_if_no_deps()
  skip_if_no_forensit()
  skip(paste("forensIT v1.2.0 perMarkerKLs() has a column-name swap",
             "in its frequency lookup; RPT collapses to 1e-20."))
  freqs <- c("1" = 0.4, "2" = 0.6)
  ped <- pedtools::nuclearPed(1)
  mp <- "3"
  # If the bug were absent, the following loop would yield exactly the
  # h1-direction oracle value computed above (chain-rule identity).
  gt <- build_geno_table(freqs)
  acc <- 0
  for (i in seq_along(gt$geno_list)) {
    for (j in seq_along(gt$geno_list)) {
      w <- gt$hwe_vec[i] * gt$hwe_vec[j]
      x <- typed_ped(ped, list("1" = gt$geno_list[[i]],
                                "2" = gt$geno_list[[j]]), freqs)
      kl <- forensIT::perMarkerKLs(x, MP = mp,
                                    frequency = list(TM = freqs))
      acc <- acc + w * as.numeric(kl$KLpedpop[[1]])
    }
  }
  mm <- mispitools::marker_model(ped, "TM", freqs)
  my_kl <- mispitools:::per_marker_kl_R(mm, poi = mp)
  expect_equal(my_kl$e_log10_lr_h1, acc, tolerance = 1e-8)
})
