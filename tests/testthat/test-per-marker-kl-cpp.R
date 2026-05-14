## F3.1 — cpp_per_marker_kl() vs per_marker_kl_R() cross-check.
##
## The C++ kernel (core::per_marker_kl) consumes the sparse joint produced
## by cpt_marker_joint and returns bidirectional KL + expected log10 LR
## with the same boundary convention as the R reference engine
## (R/r_ref_per_marker.R::per_marker_kl_R): +Inf for absolute-continuity
## violations, finite otherwise. F3.2 will wire a model-aware R wrapper;
## these tests pin the kernel itself by feeding it joint columns directly.

skip_if_no_pedtools <- function() {
  testthat::skip_if_not_installed("pedtools")
}

# Build a joint via the existing C++ wrapper, then run the F3.1 kernel
# on its (P_H1, P_H2) columns. Returns the binding's named list.
cpp_kl_from_model <- function(model, poi = NULL) {
  jt <- mispitools:::cpt_marker_joint_cpp_wrap(model, poi = poi)
  mispitools:::cpp_per_marker_kl(jt$P_H1, jt$P_H2)
}

# Compare the C++ binding's output to per_marker_kl_R(); use bit-for-bit
# tolerance for finite values and identity for ±Inf / counters.
expect_kl_equal <- function(cpp, ref, tol = 1e-12, info = NULL) {
  expect_equal(cpp$e_log10_lr_h1, ref$e_log10_lr_h1, tolerance = tol, info = info)
  expect_equal(cpp$e_log10_lr_h2, ref$e_log10_lr_h2, tolerance = tol, info = info)
  expect_equal(cpp$kl_h1_to_h2,   ref$kl_h1_to_h2,   tolerance = tol, info = info)
  expect_equal(cpp$kl_h2_to_h1,   ref$kl_h2_to_h1,   tolerance = tol, info = info)
}

test_that("cpp_per_marker_kl: identical distributions give zero KL", {
  res <- mispitools:::cpp_per_marker_kl(
    p_h1 = c(0.25, 0.25, 0.25, 0.25),
    p_h2 = c(0.25, 0.25, 0.25, 0.25)
  )
  expect_equal(res$e_log10_lr_h1, 0)
  expect_equal(res$e_log10_lr_h2, 0)
  expect_equal(res$kl_h1_to_h2, 0)
  expect_equal(res$kl_h2_to_h1, 0)
  expect_equal(res$abs_cont_violations_h1, 0L)
  expect_equal(res$abs_cont_violations_h2, 0L)
  expect_equal(res$mass_violations_h1, 0)
  expect_equal(res$mass_violations_h2, 0)
})

test_that("cpp_per_marker_kl: Bernoulli closed-form matches", {
  p <- 0.7
  q <- 0.5
  res <- mispitools:::cpp_per_marker_kl(c(p, 1 - p), c(q, 1 - q))
  e_h1 <- p * (log10(p) - log10(q)) + (1 - p) * (log10(1 - p) - log10(1 - q))
  e_h2 <- q * (log10(p) - log10(q)) + (1 - q) * (log10(1 - p) - log10(1 - q))
  expect_equal(res$e_log10_lr_h1, e_h1, tolerance = 1e-15)
  expect_equal(res$e_log10_lr_h2, e_h2, tolerance = 1e-15)
  expect_equal(res$kl_h1_to_h2, e_h1 * log(10), tolerance = 1e-15)
  expect_equal(res$kl_h2_to_h1, -e_h2 * log(10), tolerance = 1e-15)
})

test_that("cpp_per_marker_kl: H2 lacks H1 support → +Inf in H1→H2 direction", {
  res <- mispitools:::cpp_per_marker_kl(c(0.6, 0.4), c(1.0, 0.0))
  expect_identical(res$e_log10_lr_h1, Inf)
  expect_identical(res$kl_h1_to_h2, Inf)
  e_h2_want <- 1.0 * (log10(0.6) - log10(1.0))
  expect_equal(res$e_log10_lr_h2, e_h2_want, tolerance = 1e-15)
  expect_equal(res$kl_h2_to_h1, -e_h2_want * log(10), tolerance = 1e-15)
  expect_equal(res$abs_cont_violations_h2, 1L)
  expect_equal(res$abs_cont_violations_h1, 0L)
  expect_equal(res$mass_violations_h2, 0.4)
})

test_that("cpp_per_marker_kl: H1 lacks H2 support → +Inf in H2→H1 direction", {
  res <- mispitools:::cpp_per_marker_kl(c(1.0, 0.0), c(0.6, 0.4))
  expect_identical(res$e_log10_lr_h2, -Inf)
  expect_identical(res$kl_h2_to_h1, Inf)
  e_h1_want <- 1.0 * (log10(1.0) - log10(0.6))
  expect_equal(res$e_log10_lr_h1, e_h1_want, tolerance = 1e-15)
  expect_equal(res$kl_h1_to_h2, e_h1_want * log(10), tolerance = 1e-15)
  expect_equal(res$abs_cont_violations_h1, 1L)
  expect_equal(res$abs_cont_violations_h2, 0L)
  expect_equal(res$mass_violations_h1, 0.4)
})

test_that("cpp_per_marker_kl: empty joint returns zeros", {
  res <- mispitools:::cpp_per_marker_kl(numeric(0), numeric(0))
  expect_equal(res$e_log10_lr_h1, 0)
  expect_equal(res$e_log10_lr_h2, 0)
  expect_equal(res$kl_h1_to_h2, 0)
  expect_equal(res$kl_h2_to_h1, 0)
})

test_that("cpp_per_marker_kl: rejects mismatched length and negative probs", {
  expect_error(
    mispitools:::cpp_per_marker_kl(c(0.5, 0.5), c(1.0)),
    "same length"
  )
  expect_error(
    mispitools:::cpp_per_marker_kl(c(0.5, 0.5), c(-0.1, 1.1)),
    "negative probability"
  )
})

# ---------------------------------------------------------------------------
# Cross-check vs per_marker_kl_R on real joints from cpt_marker_joint_cpp.
# Mirrors the F2.2 trio K=2..4 / halfSibPed / linearPed(2) / nuclearPed(2)
# coverage so a regression in either kernel surfaces here.
# ---------------------------------------------------------------------------

test_that("cpp_per_marker_kl matches R-ref on trio K=2 mut=none", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.6)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs)
  ref <- mispitools:::per_marker_kl_R(mm)
  cpp <- cpp_kl_from_model(mm)
  # Under mutation = "none", trio K=2 is generally absolute-continuity-bounded.
  expect_kl_equal(cpp, ref, tol = 1e-12, info = "trio K=2 mut=none")
})

test_that("cpp_per_marker_kl matches R-ref on trio K=3 mut=equal R=0.005", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.3, "b" = 0.5, "c" = 0.2)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                     mutation = list(model = "equal", rate = 0.005))
  ref <- mispitools:::per_marker_kl_R(mm)
  cpp <- cpp_kl_from_model(mm)
  expect_kl_equal(cpp, ref, tol = 1e-12, info = "trio K=3 mut=equal")
})

test_that("cpp_per_marker_kl matches R-ref on halfSibPed K=3 mut=equal", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.3, "c" = 0.3)
  mm <- marker_model(pedtools::halfSibPed(), "M1", freqs,
                     mutation = list(model = "equal", rate = 0.01))
  ref <- mispitools:::per_marker_kl_R(mm)
  cpp <- cpp_kl_from_model(mm)
  expect_kl_equal(cpp, ref, tol = 1e-12, info = "halfSibPed K=3 mut=equal")
})

test_that("cpp_per_marker_kl matches R-ref on linearPed(2) mut=stepwise", {
  skip_if_no_pedtools()
  freqs <- c("12" = 0.3, "13" = 0.4, "14" = 0.3)
  mm <- marker_model(pedtools::linearPed(2), "M1", freqs,
                     mutation = list(model = "stepwise", rate = 0.005, ratio = 0.1))
  ref <- mispitools:::per_marker_kl_R(mm)
  cpp <- cpp_kl_from_model(mm)
  expect_kl_equal(cpp, ref, tol = 1e-12, info = "linearPed(2) K=3 mut=stepwise")
})

test_that("cpp_per_marker_kl matches R-ref on nuclearPed(2) mut=none (Inf branch)", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.6)
  mm <- marker_model(pedtools::nuclearPed(2), "M1", freqs)
  ref <- mispitools:::per_marker_kl_R(mm)
  cpp <- cpp_kl_from_model(mm)
  # mut=none typically yields kl_h2_to_h1 = +Inf (Mendelian-incompatible
  # states under H2 with mass under HWE marginalisation but zero under H1).
  expect_equal(cpp$kl_h1_to_h2, ref$kl_h1_to_h2, tolerance = 1e-12)
  expect_identical(cpp$kl_h2_to_h1, ref$kl_h2_to_h1)  # both +Inf
  expect_equal(cpp$e_log10_lr_h1, ref$e_log10_lr_h1, tolerance = 1e-12)
  expect_identical(cpp$e_log10_lr_h2, ref$e_log10_lr_h2)  # both -Inf
  expect_gte(cpp$abs_cont_violations_h1, 1L)
  expect_equal(cpp$abs_cont_violations_h2, 0L)
})

test_that("cpp_per_marker_kl: trio HWE K=2 p=0.5 mut=none analytic e_log10_lr_h1", {
  skip_if_no_pedtools()
  # F1.6 documents: under mut=none, p=q=0.5, trio joint gives
  # e_log10_lr_h1 = 0.625 * log10(2) bit-for-bit.
  freqs <- c("a" = 0.5, "b" = 0.5)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs)
  cpp <- cpp_kl_from_model(mm)
  expect_equal(cpp$e_log10_lr_h1, 0.625 * log10(2), tolerance = 1e-12)
})

test_that("cpp_per_marker_kl matches R-ref on trio K=4 mut=stepwise", {
  skip_if_no_pedtools()
  freqs <- c("12" = 0.1, "13" = 0.2, "14" = 0.3, "15" = 0.4)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                     mutation = list(model = "stepwise", rate = 0.005, ratio = 0.1))
  ref <- mispitools:::per_marker_kl_R(mm)
  cpp <- cpp_kl_from_model(mm)
  expect_kl_equal(cpp, ref, tol = 1e-12, info = "trio K=4 mut=stepwise")
  # All finite under stepwise R>0 (no abs-cont violations).
  expect_equal(cpp$abs_cont_violations_h1, 0L)
  expect_equal(cpp$abs_cont_violations_h2, 0L)
})
