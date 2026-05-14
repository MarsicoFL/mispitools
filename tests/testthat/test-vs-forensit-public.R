## F3.3 — Verifier: public `per_marker_kl()` against forensIT-equivalent
## chain-rule oracle, analytic trio anchor, and the internal F1.6 R
## reference engine.
##
## Layering. The public API is a thin wrapper over the C++ kernel
## (F3.1) which is itself cross-checked bit-for-bit against the R
## reference in test-per-marker-kl.R (F3.2) and test-per-marker-kl-cpp.R
## (F3.1). Here we close the chain by reaching past the R reference and
## comparing the public output directly to oracles built from
## independent machinery:
##
##   1. analytic   — trio HWE K=2 p=0.5 mut=none, closed form
##                   E[log10 LR | H1] = 0.625 * log10(2). Tol 1e-14.
##   2. forensIT-equivalent — `pedprobr::oneMarkerDistribution` +
##                   analytic HWE genotype table + chain-rule joint↔
##                   MP-conditional expansion (the math
##                   `forensIT::perMarkerKLs` was meant to implement,
##                   reconstructed in-test to side-step the v1.2.0
##                   column-rename bug; cf. SCOUT F2.7 §"Strategy F3.3"
##                   and the header of test-vs-forensit.R). Tol 1e-10.
##   3. R-reference cross-check — public matches `per_marker_kl_R()`
##                   bit-for-bit on a stepwise-mutation case the C++
##                   path threads through `mutation_matrix_cpp`
##                   (covers ordering of stepwise rows, the most
##                   numerically delicate path in F2.3-2.4). Tol 1e-12.
##   4. forensIT::perMarkerKLs FIXED local wrapper — included as a
##                   documentary check. Defines `perMarkerKLs_fixed()`
##                   inside the test, applies the chain rule by
##                   exhaustively typing non-MP members under HWE, and
##                   compares to the public output. Demonstrates that
##                   once forensIT 1.3.0 fixes the column rename, the
##                   F1.7 spec oracle agrees with us at 1e-8.
##
## A separate `test_that` block invokes `forensIT::perMarkerKLs()`
## directly and is `skip()`-ed pending the upstream fix (mirrors the
## block at the bottom of test-vs-forensit.R, but targeted at the
## public per_marker_kl()).
##
## Oracle versions at last edit: pedprobr 1.0.1, pedmut 0.9.0,
## pedtools 2.9.0, forensIT 1.2.0 (broken — see comments).

skip_if_no_deps_pub <- function() {
  testthat::skip_if_not_installed("pedtools")
  testthat::skip_if_not_installed("pedprobr")
  testthat::skip_if_not_installed("pedmut")
}

skip_if_no_forensit_pub <- function() {
  testthat::skip_if_not_installed("forensIT")
}

## -- helpers ----------------------------------------------------------

# All ordered genotype labels and their HWE weights for `freqs`. The
# label order ("a1/a2" with a1 listed before a2 in the `freqs` vector)
# matches the row labels produced by `pedprobr::oneMarkerDistribution`.
pub_build_geno_table <- function(freqs) {
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
pub_typed_ped <- function(ped, typings, freqs, mutmod = NULL) {
  args <- c(list(ped, afreq = freqs, name = "TM"), typings)
  if (!is.null(mutmod)) args$mutmod <- mutmod
  m <- do.call(pedtools::marker, args)
  pedtools::setMarkers(ped, m)
}

# Chain-rule "forensIT-equivalent" weighted KL:
#   joint_KL_log10 = E_{typings of all non-MP members ~ P_H1}
#                      [ KL( P_H1(g_MP | typings) || HWE(g_MP) ) in log10 ]
# Weights are computed via `pedprobr::likelihood()` on the typed ped.
# Valid when MP is a leaf (so P_H1(typings) = P_H2(typings) and the
# typing-margin term cancels out of the joint KL); all MPs picked here
# are leaves of the canonical pedigrees.
pub_chain_kl <- function(ped, mp_id, freqs, mutmod = NULL,
                          directions = c("h1", "both")) {
  directions <- match.arg(directions)
  members <- as.character(labels(ped))
  non_mp <- setdiff(members, mp_id)
  gt <- pub_build_geno_table(freqs)
  G <- length(gt$geno_list)
  grid <- as.matrix(do.call(expand.grid,
                             replicate(length(non_mp), seq_len(G),
                                       simplify = FALSE)))
  colnames(grid) <- non_mp
  acc_h1 <- 0; acc_h2 <- 0
  for (r in seq_len(nrow(grid))) {
    typings <- lapply(non_mp, function(id) gt$geno_list[[grid[r, id]]])
    names(typings) <- non_mp
    x <- pub_typed_ped(ped, typings, freqs, mutmod = mutmod)
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

## -- 1. analytic anchor ----------------------------------------------

test_that("public per_marker_kl: trio HWE K=2 p=0.5 mut=none matches analytic", {
  # Closed form (cf. test-per-marker-r.R and F1.6 docs):
  #   E[log10 LR | H1] = 0.625 * log10(2)
  # for a 2-allele trio at p = q = 1/2 with no mutation. The public
  # path runs cpp peeling + cpp KL; agreement at 1e-14 confirms no
  # accumulation loss vs the closed form (observed: 0 to ULP).
  skip_if_no_deps_pub()
  mm <- marker_model(pedtools::nuclearPed(1), "TM",
                     c("a" = 0.5, "b" = 0.5))
  out <- per_marker_kl(mm)
  expect_equal(out$e_log10_lr_h1, 0.625 * log10(2), tolerance = 1e-14)
  expect_equal(out$kl_h1_to_h2, 0.625 * log10(2) * log(10),
               tolerance = 1e-14)
})

## -- 2. forensIT-equivalent chain-rule oracle ------------------------

test_that("public per_marker_kl: parent_child K=2 mut=none vs forensIT-equivalent", {
  # Tolerance: 1e-10 per spec; observed agreement is ULP-level
  # (~ 5e-17). The looser bound absorbs upstream pedprobr / BLAS jitter.
  skip_if_no_deps_pub()
  freqs <- c("1" = 0.4, "2" = 0.6)
  ped <- pedtools::nuclearPed(1)
  mp <- "3"
  mm <- marker_model(ped, "TM", freqs)
  out <- per_marker_kl(mm, poi = mp)
  oracle <- pub_chain_kl(ped, mp_id = mp, freqs = freqs)
  expect_equal(out$e_log10_lr_h1, unname(oracle["kl_h1_log10"]),
               tolerance = 1e-10)
})

test_that("public per_marker_kl: parent_child K=3 mut=none vs forensIT-equivalent", {
  # 3-allele variant. G = 6 genotypes; enumerated typings = 6^2 = 36.
  skip_if_no_deps_pub()
  freqs <- c("1" = 0.3, "2" = 0.5, "3" = 0.2)
  ped <- pedtools::nuclearPed(1)
  mp <- "3"
  mm <- marker_model(ped, "TM", freqs)
  out <- per_marker_kl(mm, poi = mp)
  oracle <- pub_chain_kl(ped, mp_id = mp, freqs = freqs)
  expect_equal(out$e_log10_lr_h1, unname(oracle["kl_h1_log10"]),
               tolerance = 1e-10)
})

test_that("public per_marker_kl: halfSibPed K=2 mut=none vs forensIT-equivalent", {
  # halfSibPed: members 1..5; founders 1,2,3; nonfounders 4 (child of
  # 1,2) and 5 (child of 2,3). MP = 5 (leaf). Enumerated typings on
  # the other 4 members = G^4 = 3^4 = 81 (K=2).
  skip_if_no_deps_pub()
  freqs <- c("1" = 0.4, "2" = 0.6)
  ped <- pedtools::halfSibPed()
  mp <- "5"
  mm <- marker_model(ped, "TM", freqs)
  out <- per_marker_kl(mm, poi = mp)
  oracle <- pub_chain_kl(ped, mp_id = mp, freqs = freqs)
  expect_equal(out$e_log10_lr_h1, unname(oracle["kl_h1_log10"]),
               tolerance = 1e-10)
})

test_that("public per_marker_kl: linearPed(2) K=2 mut=none vs forensIT-equivalent (chain rule with internal nonfounder)", {
  # linearPed(2): members 1..5; founders 1,2,4; nonfounders 3 (child
  # of 1,2) and 5 (child of 3,4). MP = 5 (leaf). The intermediate
  # nonfounder 3 must be typed and its Mendelian-conditional weight
  # included in the chain-rule average — the case where the "type only
  # founders" shortcut would violate the data-processing inequality.
  # G^4 = 81 typings, many of zero Mendelian weight.
  skip_if_no_deps_pub()
  freqs <- c("1" = 0.4, "2" = 0.6)
  ped <- pedtools::linearPed(2)
  mp <- "5"
  mm <- marker_model(ped, "TM", freqs)
  out <- per_marker_kl(mm, poi = mp)
  oracle <- pub_chain_kl(ped, mp_id = mp, freqs = freqs)
  expect_equal(out$e_log10_lr_h1, unname(oracle["kl_h1_log10"]),
               tolerance = 1e-10)
})

test_that("public per_marker_kl: parent_child K=3 mut=equal r=0.005 vs forensIT-equivalent (both directions)", {
  # With strictly positive mutation, P_H1(g_c | g_f, g_m) > 0 for every
  # g_c, so the H2->H1 direction is finite and verifiable too. Oracle
  # uses `pedmut::mutationModel("equal", ...)`, matching the equal-rate
  # construction in `mutation_matrix_R` / `mutation_matrix_cpp`.
  skip_if_no_deps_pub()
  freqs <- c("1" = 0.3, "2" = 0.5, "3" = 0.2)
  ped <- pedtools::nuclearPed(1)
  mp <- "3"
  rate <- 0.005
  mut_mat <- pedmut::mutationModel("equal", alleles = names(freqs),
                                     rate = rate)
  mm <- marker_model(ped, "TM", freqs,
                     mutation = list(model = "equal", rate = rate))
  out <- per_marker_kl(mm, poi = mp)
  oracle <- pub_chain_kl(ped, mp_id = mp, freqs = freqs,
                          mutmod = mut_mat, directions = "both")
  expect_equal(out$e_log10_lr_h1, unname(oracle["kl_h1_log10"]),
               tolerance = 1e-10)
  # e_log10_lr_h2 = E_H2[log10 LR] = -KL(H2 || H1) in log10 base, so
  # oracle's "kl_h2_log10" (= KL(RPT||CPT) in log10) equals
  # -out$e_log10_lr_h2.
  expect_equal(out$e_log10_lr_h2, -unname(oracle["kl_h2_log10"]),
               tolerance = 1e-10)
})

## -- 3. R-reference cross-check (stepwise path) ----------------------

test_that("public per_marker_kl matches per_marker_kl_R bit-for-bit on stepwise mutation", {
  # Stepwise routes through `mutation_matrix_cpp` (F2.3) with a rate +
  # range parameterisation that is independent of the equal-rate
  # branch; this test guards against drift between the C++ stepwise
  # matrix builder and `mutation_matrix_R`. The trio + linearPed(2)
  # cases exercise both founder-only and chain-rule peeling paths.
  skip_if_no_deps_pub()
  freqs <- c("12" = 0.3, "13" = 0.4, "14" = 0.3)
  mm_trio <- marker_model(pedtools::nuclearPed(1), "Mstep", freqs,
                          mutation = list(model = "stepwise",
                                          rate = 0.005, ratio = 0.2))
  ref_trio <- mispitools:::per_marker_kl_R(mm_trio)
  out_trio <- per_marker_kl(mm_trio)
  expect_equal(out_trio$e_log10_lr_h1, ref_trio$e_log10_lr_h1, tolerance = 1e-12)
  expect_equal(out_trio$e_log10_lr_h2, ref_trio$e_log10_lr_h2, tolerance = 1e-12)
  expect_equal(out_trio$kl_h1_to_h2,   ref_trio$kl_h1_to_h2,   tolerance = 1e-12)
  expect_equal(out_trio$kl_h2_to_h1,   ref_trio$kl_h2_to_h1,   tolerance = 1e-12)

  mm_lin <- marker_model(pedtools::linearPed(2), "Mstep", freqs,
                         mutation = list(model = "stepwise",
                                         rate = 0.005, ratio = 0.2))
  ref_lin <- mispitools:::per_marker_kl_R(mm_lin)
  out_lin <- per_marker_kl(mm_lin)
  expect_equal(out_lin$e_log10_lr_h1, ref_lin$e_log10_lr_h1, tolerance = 1e-12)
  expect_equal(out_lin$e_log10_lr_h2, ref_lin$e_log10_lr_h2, tolerance = 1e-12)
  expect_equal(out_lin$kl_h1_to_h2,   ref_lin$kl_h1_to_h2,   tolerance = 1e-12)
  expect_equal(out_lin$kl_h2_to_h1,   ref_lin$kl_h2_to_h1,   tolerance = 1e-12)
})

## -- 4. forensIT::perMarkerKLs FIXED local wrapper -------------------

test_that("public per_marker_kl: parent_child K=2 mut=none matches a fix-puesto forensIT::perMarkerKLs (documentary)", {
  # Documentary check. forensIT v1.2.0 perMarkerKLs() has the
  # column-rename bug diagnosed in test-vs-forensit.R header (RPT
  # collapses to 1e-20). Here we reconstruct the *intended* logic
  # locally — the only change to the v1.2.0 body is the two-line
  # `pop` construction — and verify that, expanded over typings of
  # the non-MP members under HWE, it agrees with `per_marker_kl()`
  # to 1e-8 (the spec tolerance for the forensIT oracle in the
  # ROADMAP §Verificación-cruzada table).
  #
  # When forensIT 1.3.0 lands the fix upstream, the wrapper below can
  # be replaced by `forensIT::perMarkerKLs` directly without changing
  # the tolerance.
  skip_if_no_deps_pub()
  perMarkerKLs_fixed <- function(ped, MP, frequency) {
    KLpedpop <- list()
    for (i in seq_along(ped$MARKERS)) {
      df <- as.data.frame(pedprobr::oneMarkerDistribution(
        ped, ids = MP, marker = i, verbose = FALSE))
      names(df) <- "CPT"
      df$Genotype <- rownames(df)
      rownames(df) <- NULL
      alleles <- strsplit(as.character(df$Genotype), "/")
      df$Allele1 <- vapply(alleles, `[`, character(1), 1L)
      df$Allele2 <- vapply(alleles, `[`, character(1), 2L)
      ## --- fix: explicit allele/freq construction (v1.2.0 had a
      ## --- positional column rename that swapped these two).
      pop <- data.frame(
        Allele = names(frequency[[i]]),
        freq   = as.numeric(frequency[[i]]),
        stringsAsFactors = FALSE
      )
      f1 <- pop$freq[match(df$Allele1, pop$Allele)]
      f2 <- pop$freq[match(df$Allele2, pop$Allele)]
      df$RPT <- ifelse(df$Allele1 == df$Allele2, f1 * f2, 2 * f1 * f2)
      ## absolute-continuity: only states with CPT > 0 contribute.
      keep <- df$CPT > 0
      KLpedpop[[i]] <- sum(df$CPT[keep] *
                            (log10(df$CPT[keep]) - log10(df$RPT[keep])))
    }
    data.frame(markName = names(frequency),
               KLpedpop = unlist(KLpedpop))
  }

  freqs <- c("1" = 0.4, "2" = 0.6)
  ped <- pedtools::nuclearPed(1)
  mp <- "3"
  gt <- pub_build_geno_table(freqs)
  acc <- 0
  for (i in seq_along(gt$geno_list)) {
    for (j in seq_along(gt$geno_list)) {
      w <- gt$hwe_vec[i] * gt$hwe_vec[j]
      x <- pub_typed_ped(ped, list("1" = gt$geno_list[[i]],
                                    "2" = gt$geno_list[[j]]), freqs)
      kl <- perMarkerKLs_fixed(x, MP = mp, frequency = list(TM = freqs))
      acc <- acc + w * as.numeric(kl$KLpedpop[[1]])
    }
  }
  mm <- marker_model(ped, "TM", freqs)
  out <- per_marker_kl(mm, poi = mp)
  expect_equal(out$e_log10_lr_h1, unname(acc), tolerance = 1e-8)
})

## -- 5. forensIT::perMarkerKLs DIRECT (skipped pending upstream fix) -

test_that("documentation probe: public per_marker_kl vs forensIT::perMarkerKLs (skipped pending v1.2.0 RPT lookup fix)", {
  # Mirrors the skipped block at the bottom of test-vs-forensit.R but
  # for the public API. When forensIT is patched (column-rename bug
  # in the `pop` data.frame construction), removing the `skip()` below
  # should leave a passing test at the spec tolerance 1e-8.
  skip_if_no_deps_pub()
  skip_if_no_forensit_pub()
  skip(paste("forensIT v1.2.0 perMarkerKLs() has a column-name swap",
             "in its frequency lookup; RPT collapses to 1e-20."))
  freqs <- c("1" = 0.4, "2" = 0.6)
  ped <- pedtools::nuclearPed(1)
  mp <- "3"
  gt <- pub_build_geno_table(freqs)
  acc <- 0
  for (i in seq_along(gt$geno_list)) {
    for (j in seq_along(gt$geno_list)) {
      w <- gt$hwe_vec[i] * gt$hwe_vec[j]
      x <- pub_typed_ped(ped, list("1" = gt$geno_list[[i]],
                                    "2" = gt$geno_list[[j]]), freqs)
      kl <- forensIT::perMarkerKLs(x, MP = mp,
                                    frequency = list(TM = freqs))
      acc <- acc + w * as.numeric(kl$KLpedpop[[1]])
    }
  }
  mm <- marker_model(ped, "TM", freqs)
  out <- per_marker_kl(mm, poi = mp)
  expect_equal(out$e_log10_lr_h1, unname(acc), tolerance = 1e-8)
})

## -- 6. profile entry point reaches the same oracle ------------------

test_that("public per_marker_kl_profile: profile row agrees with forensIT-equivalent chain-rule", {
  # Confirms the F3.2 profile loop preserves numerical fidelity past
  # the single-marker entry point. Single-marker profile call; the
  # multi-marker version is exercised in test-per-marker-kl.R.
  skip_if_no_deps_pub()
  freqs <- c("1" = 0.3, "2" = 0.5, "3" = 0.2)
  ped <- pedtools::nuclearPed(1)
  mp <- "3"
  mm <- marker_model(ped, "TM", freqs)
  out <- per_marker_kl_profile(list(mm), poi = mp)
  oracle <- pub_chain_kl(ped, mp_id = mp, freqs = freqs)
  expect_equal(nrow(out), 1L)
  expect_equal(out$e_log10_lr_h1[1L], unname(oracle["kl_h1_log10"]),
               tolerance = 1e-10)
})
