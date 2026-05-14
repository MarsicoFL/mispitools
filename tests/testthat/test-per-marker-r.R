## F1.6 — per_marker_kl_R() / per_marker_lr_dist_R(): pure-R reference
## engine for per-marker expected log10 LR / KL bidirectional and the
## full (sparse) LR distribution. Both consume cpt_marker_joint_R from
## F1.2/F1.4/F1.5. Cross-engine oracle tests vs forensIT::perMarkerKLs
## live in F1.7.

skip_if_no_pedtools <- function() {
  testthat::skip_if_not_installed("pedtools")
}

trio_freqs <- function() c("a" = 0.4, "b" = 0.6)

trio_model <- function(freqs = trio_freqs(),
                       mutation = list(model = "none", rate = 0)) {
  skip_if_no_pedtools()
  marker_model(pedtools::nuclearPed(1), "M1", freqs, mutation = mutation)
}

# ---------------------------------------------------------------------------
# Structure / validation
# ---------------------------------------------------------------------------

test_that("per_marker_kl_R returns the documented structure", {
  mm <- trio_model()
  out <- mispitools:::per_marker_kl_R(mm)

  expect_s3_class(out, "data.frame")
  expect_named(out, c("marker", "e_log10_lr_h1", "e_log10_lr_h2",
                      "kl_h1_to_h2", "kl_h2_to_h1"))
  expect_equal(nrow(out), 1L)
  expect_identical(out$marker, "M1")
  expect_true(is.numeric(out$e_log10_lr_h1))
  expect_true(is.numeric(out$kl_h1_to_h2))
})

test_that("per_marker_lr_dist_R returns the documented structure", {
  mm <- trio_model()
  out <- mispitools:::per_marker_lr_dist_R(mm)

  expect_s3_class(out, "data.frame")
  expect_named(out, c("log10_lr", "p_h1", "p_h2"))
  expect_identical(attr(out, "marker_id"), "M1")
  expect_identical(attr(out, "poi"), "3")
  expect_true(is.numeric(out$log10_lr))
  expect_true(all(out$p_h1 >= 0))
  expect_true(all(out$p_h2 >= 0))
})

test_that("per_marker_kl_R rejects non-marker_model input", {
  expect_error(mispitools:::per_marker_kl_R(list(a = 1)),
               "must be a 'marker_model'")
})

test_that("per_marker_lr_dist_R rejects non-marker_model input", {
  expect_error(mispitools:::per_marker_lr_dist_R(list(a = 1)),
               "must be a 'marker_model'")
})

# ---------------------------------------------------------------------------
# Sums and self-consistency
# ---------------------------------------------------------------------------

test_that("p_h1 and p_h2 in LR dist each sum to 1", {
  mm <- trio_model(mutation = list(model = "equal", rate = 1e-3))
  out <- mispitools:::per_marker_lr_dist_R(mm)
  expect_equal(sum(out$p_h1), 1, tolerance = 1e-12)
  expect_equal(sum(out$p_h2), 1, tolerance = 1e-12)
})

test_that("p_h1 and p_h2 sum to 1 also under mutation = none", {
  mm <- trio_model()
  out <- mispitools:::per_marker_lr_dist_R(mm)
  expect_equal(sum(out$p_h1), 1, tolerance = 1e-12)
  expect_equal(sum(out$p_h2), 1, tolerance = 1e-12)
})

test_that("kl_h1_to_h2 = e_log10_lr_h1 * ln(10)  (self-consistency)", {
  mm <- trio_model(mutation = list(model = "equal", rate = 1e-3))
  out <- mispitools:::per_marker_kl_R(mm)
  expect_equal(out$kl_h1_to_h2, out$e_log10_lr_h1 * log(10),
               tolerance = 1e-12)
})

test_that("kl_h2_to_h1 = -e_log10_lr_h2 * ln(10)  (self-consistency)", {
  mm <- trio_model(mutation = list(model = "equal", rate = 1e-3))
  out <- mispitools:::per_marker_kl_R(mm)
  expect_equal(out$kl_h2_to_h1, -out$e_log10_lr_h2 * log(10),
               tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Boundary: mutation = none yields divergent KL(H2 || H1)
# ---------------------------------------------------------------------------

test_that("mutation = none: kl_h2_to_h1 = +Inf and e_log10_lr_h2 = -Inf", {
  out <- mispitools:::per_marker_kl_R(trio_model())
  expect_identical(out$kl_h2_to_h1, Inf)
  expect_identical(out$e_log10_lr_h2, -Inf)
})

test_that("mutation = none: kl_h1_to_h2 stays finite and non-negative", {
  out <- mispitools:::per_marker_kl_R(trio_model())
  expect_true(is.finite(out$kl_h1_to_h2))
  expect_gte(out$kl_h1_to_h2, 0)
  expect_true(is.finite(out$e_log10_lr_h1))
  expect_gte(out$e_log10_lr_h1, 0)
})

test_that("mutation = none: LR dist has rows with log10_lr = -Inf", {
  out <- mispitools:::per_marker_lr_dist_R(trio_model())
  expect_true(any(is.infinite(out$log10_lr) & out$log10_lr < 0))
  ## All -Inf rows must have p_h1 = 0 and p_h2 > 0.
  neg_inf <- is.infinite(out$log10_lr) & out$log10_lr < 0
  expect_true(all(out$p_h1[neg_inf] == 0))
  expect_true(all(out$p_h2[neg_inf] > 0))
})

# ---------------------------------------------------------------------------
# Boundary: mutation = equal > 0 → both KLs finite, both non-negative
# ---------------------------------------------------------------------------

test_that("mutation = equal: both KLs finite and non-negative", {
  out <- mispitools:::per_marker_kl_R(
    trio_model(mutation = list(model = "equal", rate = 1e-3))
  )
  expect_true(is.finite(out$kl_h1_to_h2))
  expect_true(is.finite(out$kl_h2_to_h1))
  expect_gte(out$kl_h1_to_h2, 0)
  expect_gte(out$kl_h2_to_h1, 0)
  expect_true(is.finite(out$e_log10_lr_h1))
  expect_true(is.finite(out$e_log10_lr_h2))
  expect_gte(out$e_log10_lr_h1, 0)
  expect_lte(out$e_log10_lr_h2, 0)
})

test_that("mutation = stepwise: both KLs finite and non-negative", {
  mm <- marker_model(
    pedtools::nuclearPed(1), "M1",
    freqs = c("10" = 0.2, "11" = 0.3, "12" = 0.5),
    mutation = list(model = "stepwise", rate = 1e-3, ratio = 0.1)
  )
  out <- mispitools:::per_marker_kl_R(mm)
  expect_true(is.finite(out$kl_h1_to_h2))
  expect_true(is.finite(out$kl_h2_to_h1))
  expect_gte(out$kl_h1_to_h2, 0)
  expect_gte(out$kl_h2_to_h1, 0)
})

# ---------------------------------------------------------------------------
# Analytic: parent-child, mutation = none, uniform 2-allele freqs.
# Derivation (see r_ref_per_marker.R header):
#   KL(H1 || H2) = E_{g1,g2 ~ HWE}[ KL(P(.|g1,g2) || P_HWE(.)) ]
# For p = q = 0.5 this evaluates to 0.625 * log(2), so
#   e_log10_lr_h1 = 0.625 * log10(2)
# ---------------------------------------------------------------------------

test_that("analytic trio p=q=0.5 no-mutation: e_log10_lr_h1 = 0.625 * log10(2)", {
  out <- mispitools:::per_marker_kl_R(
    trio_model(freqs = c("a" = 0.5, "b" = 0.5))
  )
  expect_equal(out$e_log10_lr_h1, 0.625 * log10(2), tolerance = 1e-12)
  expect_equal(out$kl_h1_to_h2, 0.625 * log(2), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Continuity: as mutation rate -> 0 with `equal`, kl_h2_to_h1 diverges
# (the regularizer that protects against -Inf states shrinks).
# ---------------------------------------------------------------------------

test_that("kl_h2_to_h1 grows as equal-rate mutation rate -> 0", {
  freqs <- c("a" = 0.4, "b" = 0.6)
  kl_at <- function(rate) {
    mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                       mutation = list(model = "equal", rate = rate))
    mispitools:::per_marker_kl_R(mm)$kl_h2_to_h1
  }
  expect_gt(kl_at(1e-5), kl_at(1e-3))
  expect_gt(kl_at(1e-3), kl_at(1e-1))
})

# ---------------------------------------------------------------------------
# Aggregation: rows with the same log10_lr collapse; sums preserved.
# ---------------------------------------------------------------------------

test_that("aggregate = TRUE never increases row count vs aggregate = FALSE", {
  mm <- trio_model(mutation = list(model = "equal", rate = 1e-3))
  raw <- mispitools:::per_marker_lr_dist_R(mm, aggregate = FALSE)
  agg <- mispitools:::per_marker_lr_dist_R(mm, aggregate = TRUE)
  expect_lte(nrow(agg), nrow(raw))
  expect_equal(sum(agg$p_h1), sum(raw$p_h1), tolerance = 1e-12)
  expect_equal(sum(agg$p_h2), sum(raw$p_h2), tolerance = 1e-12)
})

test_that("aggregated log10_lr keys are unique and sorted ascending", {
  mm <- trio_model(mutation = list(model = "equal", rate = 1e-3))
  agg <- mispitools:::per_marker_lr_dist_R(mm)
  expect_false(any(duplicated(agg$log10_lr)))
  expect_identical(agg$log10_lr, sort(agg$log10_lr))
})

# ---------------------------------------------------------------------------
# Consistency with the joint table the engine consumes.
# ---------------------------------------------------------------------------

test_that("LR dist log10_lr matches log10(P_H1 / P_H2) of joint table rows", {
  mm <- trio_model(mutation = list(model = "equal", rate = 1e-3))
  joint <- mispitools:::cpt_marker_joint_R(mm)
  raw <- mispitools:::per_marker_lr_dist_R(mm, aggregate = FALSE)
  expected <- log10(joint$P_H1) - log10(joint$P_H2)
  expect_equal(raw$log10_lr, expected, tolerance = 1e-12)
})

test_that("e_log10_lr_h1 equals manually-computed E_H1[log10 LR]", {
  mm <- trio_model(mutation = list(model = "equal", rate = 1e-3))
  joint <- mispitools:::cpt_marker_joint_R(mm)
  ll <- log10(joint$P_H1) - log10(joint$P_H2)
  manual <- sum(joint$P_H1[joint$P_H1 > 0] * ll[joint$P_H1 > 0])
  out <- mispitools:::per_marker_kl_R(mm)
  expect_equal(out$e_log10_lr_h1, manual, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# POI propagation
# ---------------------------------------------------------------------------

test_that("explicit `poi` flows through per_marker_lr_dist_R", {
  mm <- trio_model()
  out <- mispitools:::per_marker_lr_dist_R(mm, poi = "3")
  expect_identical(attr(out, "poi"), "3")
})

test_that("changing POI changes the LR distribution (parent vs child)", {
  mm <- trio_model(mutation = list(model = "equal", rate = 1e-3))
  out_child <- mispitools:::per_marker_kl_R(mm, poi = "3")
  out_father <- mispitools:::per_marker_kl_R(mm, poi = "1")
  ## With symmetric topology these should match numerically when freqs
  ## are uniform — but freqs = (0.4, 0.6) is asymmetric only via labels,
  ## not topology. Father and child play different structural roles, so
  ## the KLs should not coincide bit-for-bit.
  expect_false(isTRUE(all.equal(out_child$kl_h1_to_h2,
                                out_father$kl_h1_to_h2,
                                tolerance = 1e-8)))
})
