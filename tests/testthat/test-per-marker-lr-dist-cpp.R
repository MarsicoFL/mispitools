## F4.1 — cpp_per_marker_lr_dist() vs per_marker_lr_dist_R() cross-check.
##
## The C++ kernel (core::per_marker_lr_dist) consumes the sparse joint
## produced by cpt_marker_joint and returns the sparse-sorted
## (log10_lr, p_h1, p_h2) atoms with the same boundary convention and
## aggregation as the R reference engine
## (R/r_ref_per_marker.R::per_marker_lr_dist_R). F4.4 will wire a
## model-aware R wrapper (lr_distribution()); these tests pin the kernel
## by feeding it joint columns directly, mirroring the F3.1 KL chain.

skip_if_no_pedtools <- function() {
  testthat::skip_if_not_installed("pedtools")
}

# Build a joint via the existing C++ wrapper, then run the F4.1 kernel
# on its (P_H1, P_H2) columns.
cpp_lrdist_from_model <- function(model, poi = NULL, aggregate = TRUE) {
  jt <- mispitools:::cpt_marker_joint_cpp_wrap(model, poi = poi)
  mispitools:::cpp_per_marker_lr_dist(jt$P_H1, jt$P_H2, aggregate = aggregate)
}

# Compare an aggregated cpp distribution to per_marker_lr_dist_R(); both
# are sorted ascending by log10_lr, so a parallel scan is well-defined.
# all.equal (used by expect_equal) treats Inf == Inf correctly.
expect_lrdist_equal <- function(cpp, ref, tol = 1e-12, info = NULL) {
  expect_equal(length(cpp$log10_lr), nrow(ref), info = info)
  expect_equal(cpp$log10_lr, ref$log10_lr, tolerance = tol, info = info)
  expect_equal(cpp$p_h1, ref$p_h1, tolerance = tol, info = info)
  expect_equal(cpp$p_h2, ref$p_h2, tolerance = tol, info = info)
}

# ---------------------------------------------------------------------------
# Direct kernel behaviour on synthetic joints.
# ---------------------------------------------------------------------------

test_that("cpp_per_marker_lr_dist: Bernoulli, sorted ascending, sums to 1", {
  res <- mispitools:::cpp_per_marker_lr_dist(c(0.7, 0.3), c(0.5, 0.5))
  expect_equal(length(res$log10_lr), 2L)
  expect_true(all(diff(res$log10_lr) > 0))
  # ascending: the (p1=0.3) atom first, (p1=0.7) atom second.
  expect_equal(res$log10_lr[1], log10(0.3) - log10(0.5), tolerance = 1e-15)
  expect_equal(res$log10_lr[2], log10(0.7) - log10(0.5), tolerance = 1e-15)
  expect_equal(sum(res$p_h1), 1, tolerance = 1e-15)
  expect_equal(sum(res$p_h2), 1, tolerance = 1e-15)
  expect_false(res$has_pos_inf)
  expect_false(res$has_neg_inf)
})

test_that("cpp_per_marker_lr_dist: equal-log10_lr atoms collapse", {
  res <- mispitools:::cpp_per_marker_lr_dist(c(0.2, 0.2, 0.6),
                                             c(0.4, 0.4, 0.2))
  expect_equal(length(res$log10_lr), 2L)
  expect_equal(res$log10_lr[1], log10(0.2) - log10(0.4), tolerance = 1e-15)
  expect_equal(res$p_h1, c(0.4, 0.6), tolerance = 1e-15)
  expect_equal(res$p_h2, c(0.8, 0.2), tolerance = 1e-15)
})

test_that("cpp_per_marker_lr_dist: +Inf atom sorts last, flag set", {
  res <- mispitools:::cpp_per_marker_lr_dist(c(0.6, 0.4), c(1.0, 0.0))
  expect_equal(length(res$log10_lr), 2L)
  expect_identical(res$log10_lr[2], Inf)
  expect_equal(res$p_h1[2], 0.4, tolerance = 1e-15)
  expect_equal(res$p_h2[2], 0.0)
  expect_true(res$has_pos_inf)
  expect_false(res$has_neg_inf)
})

test_that("cpp_per_marker_lr_dist: -Inf atom sorts first, flag set", {
  res <- mispitools:::cpp_per_marker_lr_dist(c(0.0, 0.6), c(0.4, 0.6))
  expect_identical(res$log10_lr[1], -Inf)
  expect_true(res$has_neg_inf)
  expect_false(res$has_pos_inf)
})

test_that("cpp_per_marker_lr_dist: empty joint returns empty", {
  res <- mispitools:::cpp_per_marker_lr_dist(numeric(0), numeric(0))
  expect_equal(length(res$log10_lr), 0L)
  expect_false(res$has_pos_inf)
  expect_false(res$has_neg_inf)
})

test_that("cpp_per_marker_lr_dist: rejects mismatch and negative probs", {
  expect_error(
    mispitools:::cpp_per_marker_lr_dist(c(0.5, 0.5), c(1.0)),
    "same length"
  )
  expect_error(
    mispitools:::cpp_per_marker_lr_dist(c(0.5, 0.5), c(-0.1, 1.1)),
    "negative probability"
  )
})

test_that("cpp_per_marker_lr_dist: aggregate=FALSE keeps raw atoms, sorting recovers aggregate", {
  p1 <- c(0.2, 0.6, 0.2)
  p2 <- c(0.4, 0.2, 0.4)
  raw <- mispitools:::cpp_per_marker_lr_dist(p1, p2, aggregate = FALSE)
  agg <- mispitools:::cpp_per_marker_lr_dist(p1, p2, aggregate = TRUE)
  expect_equal(length(raw$log10_lr), 3L)   # no collapse
  expect_equal(length(agg$log10_lr), 2L)   # two distinct keys
  # Manually aggregate the raw atoms → must equal the aggregated path.
  d <- data.frame(k = raw$log10_lr, p1 = raw$p_h1, p2 = raw$p_h2)
  d <- d[order(d$k), , drop = FALSE]
  s1 <- tapply(d$p1, d$k, sum)
  s2 <- tapply(d$p2, d$k, sum)
  expect_equal(unname(as.numeric(s1)), agg$p_h1, tolerance = 1e-15)
  expect_equal(unname(as.numeric(s2)), agg$p_h2, tolerance = 1e-15)
})

# ---------------------------------------------------------------------------
# Cross-check vs per_marker_lr_dist_R on real joints from cpt_marker_joint.
# Mirrors the F3.1 trio K=2..4 / halfSibPed / linearPed(2) / nuclearPed(2)
# coverage so a regression in either kernel surfaces here.
# ---------------------------------------------------------------------------

test_that("cpp_per_marker_lr_dist matches R-ref on trio K=2 mut=none", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.6)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs)
  ref <- mispitools:::per_marker_lr_dist_R(mm)
  cpp <- cpp_lrdist_from_model(mm)
  expect_lrdist_equal(cpp, ref, info = "trio K=2 mut=none")
})

test_that("cpp_per_marker_lr_dist matches R-ref on trio K=3 mut=equal", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.3, "b" = 0.5, "c" = 0.2)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                     mutation = list(model = "equal", rate = 0.005))
  ref <- mispitools:::per_marker_lr_dist_R(mm)
  cpp <- cpp_lrdist_from_model(mm)
  expect_lrdist_equal(cpp, ref, info = "trio K=3 mut=equal")
})

test_that("cpp_per_marker_lr_dist matches R-ref on halfSibPed K=3 mut=equal", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.3, "c" = 0.3)
  mm <- marker_model(pedtools::halfSibPed(), "M1", freqs,
                     mutation = list(model = "equal", rate = 0.01))
  ref <- mispitools:::per_marker_lr_dist_R(mm)
  cpp <- cpp_lrdist_from_model(mm)
  expect_lrdist_equal(cpp, ref, info = "halfSibPed K=3 mut=equal")
})

test_that("cpp_per_marker_lr_dist matches R-ref on linearPed(2) mut=stepwise", {
  skip_if_no_pedtools()
  freqs <- c("12" = 0.3, "13" = 0.4, "14" = 0.3)
  mm <- marker_model(pedtools::linearPed(2), "M1", freqs,
                     mutation = list(model = "stepwise", rate = 0.005, ratio = 0.1))
  ref <- mispitools:::per_marker_lr_dist_R(mm)
  cpp <- cpp_lrdist_from_model(mm)
  expect_lrdist_equal(cpp, ref, info = "linearPed(2) K=3 mut=stepwise")
})

test_that("cpp_per_marker_lr_dist matches R-ref on nuclearPed(2) mut=none (Inf atoms)", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.6)
  mm <- marker_model(pedtools::nuclearPed(2), "M1", freqs)
  ref <- mispitools:::per_marker_lr_dist_R(mm)
  cpp <- cpp_lrdist_from_model(mm)
  expect_lrdist_equal(cpp, ref, info = "nuclearPed(2) mut=none")
  # mut=none on this pedigree yields -Inf atoms (H1 lacks Mendelian-
  # incompatible states that H2 supports under HWE marginalisation).
  expect_true(cpp$has_neg_inf)
  expect_identical(cpp$log10_lr[1], -Inf)
})

test_that("cpp_per_marker_lr_dist matches R-ref on trio K=4 mut=stepwise", {
  skip_if_no_pedtools()
  freqs <- c("12" = 0.1, "13" = 0.2, "14" = 0.3, "15" = 0.4)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                     mutation = list(model = "stepwise", rate = 0.005, ratio = 0.1))
  ref <- mispitools:::per_marker_lr_dist_R(mm)
  cpp <- cpp_lrdist_from_model(mm)
  expect_lrdist_equal(cpp, ref, info = "trio K=4 mut=stepwise")
  expect_false(cpp$has_pos_inf)
  expect_false(cpp$has_neg_inf)
})

test_that("cpp_per_marker_lr_dist: sum of p_h1 / p_h2 is 1 on a real joint", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.3, "b" = 0.5, "c" = 0.2)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                     mutation = list(model = "equal", rate = 0.005))
  cpp <- cpp_lrdist_from_model(mm)
  expect_equal(sum(cpp$p_h1), 1, tolerance = 1e-12)
  expect_equal(sum(cpp$p_h2), 1, tolerance = 1e-12)
})

test_that("cpp_per_marker_lr_dist: mean log10 LR consistent with per_marker_kl", {
  skip_if_no_pedtools()
  # Self-consistency anchor (SCOUT_DNAtools F4.2 invariant precursor):
  # E[log10 LR | H1] = sum p_h1 * log10_lr must equal per_marker_kl's
  # e_log10_lr_h1 (finite case, mut=equal so no Inf atoms).
  freqs <- c("a" = 0.3, "b" = 0.5, "c" = 0.2)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                     mutation = list(model = "equal", rate = 0.005))
  cpp <- cpp_lrdist_from_model(mm)
  kl  <- mispitools:::per_marker_kl_R(mm)
  expect_equal(sum(cpp$p_h1 * cpp$log10_lr), kl$e_log10_lr_h1,
               tolerance = 1e-12)
  expect_equal(sum(cpp$p_h2 * cpp$log10_lr), kl$e_log10_lr_h2,
               tolerance = 1e-12)
})
