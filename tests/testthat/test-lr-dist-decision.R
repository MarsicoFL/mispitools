## F4.3 — decision-theoretic primitives (core) + S3 layer.
##
## The C++ kernel (core::lr_dist_summary / lr_dist_quantile /
## decision_rates / roc_curve / choose_threshold_weighted) is pinned
## against analytic closed forms and the R reference
## (R/r_ref_per_marker.R::*_R) on both synthetic distributions and ones
## derived from a real marker_model via the F4.1 chain
## (cpt_marker_joint_cpp_wrap -> cpp_per_marker_lr_dist). The R S3
## methods (as_lr_dist / summary / plot / print) are smoke-tested.

skip_if_no_pedtools <- function() {
  testthat::skip_if_not_installed("pedtools")
}

# Reference distribution with known closed forms.
ref_d <- list(
  log10_lr = c(-1, 0, 2),
  p_h1 = c(0.1, 0.3, 0.6),
  p_h2 = c(0.6, 0.3, 0.1)
)

cpp_lrdist_from_model <- function(model, poi = NULL) {
  jt <- mispitools:::cpt_marker_joint_cpp_wrap(model, poi = poi)
  mispitools:::cpp_per_marker_lr_dist(jt$P_H1, jt$P_H2, aggregate = TRUE)
}

# ---------------------------------------------------------------------------
# Closed-form anchors.
# ---------------------------------------------------------------------------

test_that("lr_dist_summary: analytic moments", {
  s <- mispitools:::cpp_lr_dist_summary(
    ref_d$log10_lr, ref_d$p_h1, ref_d$p_h2)
  expect_equal(s$mean_h1, 1.1, tolerance = 1e-12)
  expect_equal(s$mean_h2, -0.4, tolerance = 1e-12)
  expect_equal(s$var_h1, 1.29, tolerance = 1e-12)
  expect_equal(s$mass_h1, 1, tolerance = 1e-12)
  expect_equal(s$mass_h2, 1, tolerance = 1e-12)
  expect_false(s$has_pos_inf)
  expect_false(s$has_neg_inf)
})

test_that("decision_rates: == threshold atom is indeterminate", {
  v <- mispitools:::cpp_lr_dist_decision_rates(
    ref_d$log10_lr, ref_d$p_h1, ref_d$p_h2, threshold = 0)
  expect_equal(v$fnr, 0.1, tolerance = 1e-12)
  expect_equal(v$fpr, 0.1, tolerance = 1e-12)
  expect_equal(v$tpr, 0.6, tolerance = 1e-12)
  expect_equal(v$tnr, 0.6, tolerance = 1e-12)
  # tpr != 1 - fnr because the log10_lr == 0 atom is excluded.
  expect_false(isTRUE(all.equal(v$tpr, 1 - v$fnr)))
})

test_that("roc_curve: concordance AUC closed form", {
  r <- mispitools:::cpp_lr_dist_roc(
    ref_d$log10_lr, ref_d$p_h1, ref_d$p_h2)
  expect_equal(r$auc, 0.825, tolerance = 1e-12)
  expect_equal(length(r$threshold), 3L)
  # AUC equals trapezoid area under the (fpr, tpr) step curve.
})

test_that("lr_dist_quantile: type-1 inverse CDF", {
  q <- mispitools:::cpp_lr_dist_quantile(
    ref_d$log10_lr, ref_d$p_h1, ref_d$p_h2,
    probs = c(0.1, 0.4, 0.5), under_h1 = TRUE)
  expect_equal(q, c(-1, 0, 2), tolerance = 1e-12)
})

test_that("choose_threshold_weighted: weighted-Euclidean optimum", {
  v <- mispitools:::cpp_lr_dist_choose_threshold(
    ref_d$log10_lr, ref_d$p_h1, ref_d$p_h2, weight = 1)
  expect_equal(v$threshold, 0, tolerance = 1e-12)
  expect_equal(v$distance, sqrt(0.02), tolerance = 1e-12)
  expect_equal(v$fpr, 0.1, tolerance = 1e-12)
  expect_equal(v$fnr, 0.1, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# C++ kernel bit-for-bit vs R reference (synthetic).
# ---------------------------------------------------------------------------

test_that("cpp == R-ref on a synthetic distribution", {
  set.seed(20260516)
  lr <- sort(c(-Inf, runif(8, -3, 3)))
  p1 <- runif(9); p1 <- p1 / sum(p1)
  p2 <- runif(9); p2 <- p2 / sum(p2)
  # -Inf atom: P_H1 = 0 there (boundary convention).
  p1[1] <- 0; p1 <- p1 / sum(p1)
  d <- list(log10_lr = lr, p_h1 = p1, p_h2 = p2)

  s_cpp <- mispitools:::cpp_lr_dist_summary(lr, p1, p2)
  s_ref <- mispitools:::lr_dist_summary_R(d)
  expect_equal(s_cpp$mean_h1, s_ref$mean_h1, tolerance = 1e-12)
  expect_equal(s_cpp$mean_h2, s_ref$mean_h2, tolerance = 1e-12)
  expect_equal(s_cpp$var_h1, s_ref$var_h1, tolerance = 1e-12)
  expect_equal(s_cpp$var_h2, s_ref$var_h2, tolerance = 1e-12)
  expect_equal(s_cpp$has_neg_inf, s_ref$has_neg_inf)

  probs <- c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)
  expect_equal(
    mispitools:::cpp_lr_dist_quantile(lr, p1, p2, probs, TRUE),
    mispitools:::lr_dist_quantile_R(d, probs, TRUE),
    tolerance = 1e-12)
  expect_equal(
    mispitools:::cpp_lr_dist_quantile(lr, p1, p2, probs, FALSE),
    mispitools:::lr_dist_quantile_R(d, probs, FALSE),
    tolerance = 1e-12)

  for (t in c(-1.5, 0, 1.2)) {
    rc <- mispitools:::cpp_lr_dist_decision_rates(lr, p1, p2, t)
    rr <- mispitools:::decision_rates_R(d, t)
    expect_equal(rc$fpr, rr$fpr, tolerance = 1e-12)
    expect_equal(rc$fnr, rr$fnr, tolerance = 1e-12)
    expect_equal(rc$mcc, rr$mcc, tolerance = 1e-12)
  }

  roc_cpp <- mispitools:::cpp_lr_dist_roc(lr, p1, p2)
  roc_ref <- mispitools:::roc_curve_R(d)
  expect_equal(roc_cpp$auc, roc_ref$auc, tolerance = 1e-12)
  expect_equal(roc_cpp$fpr, roc_ref$fpr, tolerance = 1e-12)
  expect_equal(roc_cpp$tpr, roc_ref$tpr, tolerance = 1e-12)

  ct_cpp <- mispitools:::cpp_lr_dist_choose_threshold(lr, p1, p2, 10)
  ct_ref <- mispitools:::choose_threshold_weighted_R(d, 10)
  expect_equal(ct_cpp$threshold, ct_ref$threshold, tolerance = 1e-12)
  expect_equal(ct_cpp$distance, ct_ref$distance, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# C++ kernel bit-for-bit vs R reference (model-derived).
# ---------------------------------------------------------------------------

test_that("cpp == R-ref on a model-derived LR distribution", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  freqs <- c(A = 0.2, B = 0.3, C = 0.5)
  for (mut in list(list(model = "none", rate = 0),
                   list(model = "equal", rate = 0.005))) {
    model <- marker_model(ped, "M1", freqs, mutation = mut)
    d <- cpp_lrdist_from_model(model)
    dd <- list(log10_lr = d$log10_lr, p_h1 = d$p_h1, p_h2 = d$p_h2)

    s_cpp <- mispitools:::cpp_lr_dist_summary(d$log10_lr, d$p_h1, d$p_h2)
    s_ref <- mispitools:::lr_dist_summary_R(dd)
    expect_equal(s_cpp$mean_h1, s_ref$mean_h1, tolerance = 1e-12,
                 info = mut$model)
    expect_equal(s_cpp$var_h2, s_ref$var_h2, tolerance = 1e-12,
                 info = mut$model)

    # E[log10 LR | Hk] must agree with the F3 KL engine.
    kl <- mispitools:::cpp_per_marker_kl(d$p_h1, d$p_h2)
    expect_equal(s_cpp$mean_h1, kl$e_log10_lr_h1, tolerance = 1e-10,
                 info = mut$model)
    expect_equal(s_cpp$mean_h2, kl$e_log10_lr_h2, tolerance = 1e-10,
                 info = mut$model)

    roc_cpp <- mispitools:::cpp_lr_dist_roc(d$log10_lr, d$p_h1, d$p_h2)
    roc_ref <- mispitools:::roc_curve_R(dd)
    expect_equal(roc_cpp$auc, roc_ref$auc, tolerance = 1e-12,
                 info = mut$model)
    # AUC is a probability.
    expect_gte(roc_cpp$auc, 0)
    expect_lte(roc_cpp$auc, 1)
  }
})

# ---------------------------------------------------------------------------
# Error paths.
# ---------------------------------------------------------------------------

test_that("kernels reject malformed input", {
  expect_error(
    mispitools:::cpp_lr_dist_summary(c(0, 1), 0.5, c(0.5, 0.5)),
    "same length")
  expect_error(
    mispitools:::cpp_lr_dist_summary(c(0, 1), c(-0.1, 1), c(0.5, 0.5)),
    "non-negative")
  expect_error(
    mispitools:::cpp_lr_dist_choose_threshold(
      ref_d$log10_lr, ref_d$p_h1, ref_d$p_h2, weight = 0),
    "positive")
  expect_error(
    mispitools:::cpp_lr_dist_quantile(
      c(0, 1), c(0, 0), c(0.5, 0.5), probs = 0.5, under_h1 = TRUE),
    "no active support")
})

# ---------------------------------------------------------------------------
# S3 layer.
# ---------------------------------------------------------------------------

test_that("as_lr_dist constructs and validates", {
  x <- as_lr_dist(as.data.frame(ref_d))
  expect_s3_class(x, "lr_dist")
  expect_identical(as_lr_dist(x), x)
  expect_error(as_lr_dist(list(log10_lr = 1)), "components")
  expect_error(
    as_lr_dist(list(log10_lr = c(1, 2), p_h1 = 1, p_h2 = c(1, 0))),
    "same length")
  expect_error(
    as_lr_dist(list(log10_lr = 1, p_h1 = -1, p_h2 = 1)),
    "non-negative")
})

test_that("summary.lr_dist returns the decision summary", {
  s <- summary(as_lr_dist(as.data.frame(ref_d)))
  expect_s3_class(s, "summary.lr_dist")
  expect_equal(s$mean_h1, 1.1, tolerance = 1e-12)
  expect_equal(unname(s$quantiles_h1["50%"]), 2, tolerance = 1e-12)
  expect_equal(s$auc, 0.825, tolerance = 1e-12)
  expect_output(print(s), "AUC")
})

test_that("plot.lr_dist returns a ggplot and drops Inf atoms", {
  skip_if_not_installed("ggplot2")
  g <- plot(as_lr_dist(as.data.frame(ref_d)))
  expect_s3_class(g, "ggplot")
  inf_d <- as_lr_dist(list(log10_lr = c(-Inf, 0, 1),
                           p_h1 = c(0, 0.5, 0.5),
                           p_h2 = c(0.5, 0.25, 0.25)))
  expect_message(plot(inf_d), "infinite atom")
})

test_that("print.lr_dist is non-destructive", {
  x <- as_lr_dist(as.data.frame(ref_d))
  expect_output(print(x), "lr_dist")
  expect_identical(print(x), x)
})
