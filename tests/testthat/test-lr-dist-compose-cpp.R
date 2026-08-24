## F4.2 — cpp_lr_dist_compose() vs lr_dist_compose_R() cross-check.
##
## core::lr_dist_compose convolves independent per-feature LR
## distributions (the total log10 LR is the sum of per-feature log10
## LRs under conditional independence). Exact mode reproduces the R
## reference (R/r_ref_per_marker.R::lr_dist_compose_R) bit-for-bit to
## 1e-12; the grid heuristic preserves total mass and the composed mean
## exactly with O(delta) shape error. The model-aware R wrapper
## (lr_distribution()) arrives in F4.4; these tests pin the kernel by
## feeding it per-feature distribution lists directly.

skip_if_no_pedtools <- function() {
  testthat::skip_if_not_installed("pedtools")
}

# Per-marker LR distribution for a model as a plain list(log10_lr,p_h1,p_h2),
# the input shape cpp_lr_dist_compose expects.
lrdist_of <- function(model, poi = NULL) {
  jt <- mispitools:::cpt_marker_joint_cpp_wrap(model, poi = poi)
  d <- mispitools:::cpp_per_marker_lr_dist(jt$P_H1, jt$P_H2, aggregate = TRUE)
  list(log10_lr = d$log10_lr, p_h1 = d$p_h1, p_h2 = d$p_h2)
}

as_df <- function(d) {
  data.frame(log10_lr = d$log10_lr, p_h1 = d$p_h1, p_h2 = d$p_h2)
}

expect_dist_equal <- function(cpp, ref, tol = 1e-12, info = NULL) {
  expect_equal(length(cpp$log10_lr), nrow(ref), info = info)
  expect_equal(cpp$log10_lr, ref$log10_lr, tolerance = tol, info = info)
  expect_equal(cpp$p_h1, ref$p_h1, tolerance = tol, info = info)
  expect_equal(cpp$p_h2, ref$p_h2, tolerance = tol, info = info)
}

# ---------------------------------------------------------------------------
# Direct kernel behaviour on synthetic distributions.
# ---------------------------------------------------------------------------

test_that("exact compose of two Bernoulli features = closed-form convolution", {
  d1 <- list(log10_lr = c(-1, 1), p_h1 = c(0.3, 0.7), p_h2 = c(0.6, 0.4))
  d2 <- list(log10_lr = c(0, 2),  p_h1 = c(0.5, 0.5), p_h2 = c(0.2, 0.8))
  res <- mispitools:::cpp_lr_dist_compose(list(d1, d2))
  # support = {-1,1,3} (1 appears twice: -1+2 and 1+0 → collapse)
  expect_equal(res$log10_lr, c(-1, 1, 3), tolerance = 1e-12)
  expect_equal(res$p_h1,
               c(0.3 * 0.5, 0.3 * 0.5 + 0.7 * 0.5, 0.7 * 0.5),
               tolerance = 1e-12)
  expect_equal(res$p_h2,
               c(0.6 * 0.2, 0.6 * 0.8 + 0.4 * 0.2, 0.4 * 0.8),
               tolerance = 1e-12)
  expect_equal(sum(res$p_h1), 1, tolerance = 1e-12)
  expect_equal(sum(res$p_h2), 1, tolerance = 1e-12)
})

test_that("empty list composes to the identity delta", {
  res <- mispitools:::cpp_lr_dist_compose(list())
  expect_equal(res$log10_lr, 0)
  expect_equal(res$p_h1, 1)
  expect_equal(res$p_h2, 1)
})

test_that("idempotence: compose([d]) == d", {
  d <- list(log10_lr = c(-0.7, 0.2, 1.4),
            p_h1 = c(0.1, 0.5, 0.4),
            p_h2 = c(0.5, 0.3, 0.2))
  res <- mispitools:::cpp_lr_dist_compose(list(d))
  expect_equal(res$log10_lr, d$log10_lr, tolerance = 1e-15)
  expect_equal(res$p_h1, d$p_h1, tolerance = 1e-15)
  expect_equal(res$p_h2, d$p_h2, tolerance = 1e-15)
})

test_that("commutativity: compose([d1,d2]) == compose([d2,d1])", {
  d1 <- list(log10_lr = c(-1.3, 0.4, 2.1),
             p_h1 = c(0.2, 0.5, 0.3), p_h2 = c(0.4, 0.4, 0.2))
  d2 <- list(log10_lr = c(0.1, 0.9),
             p_h1 = c(0.55, 0.45), p_h2 = c(0.25, 0.75))
  a <- mispitools:::cpp_lr_dist_compose(list(d1, d2))
  b <- mispitools:::cpp_lr_dist_compose(list(d2, d1))
  expect_equal(a$log10_lr, b$log10_lr, tolerance = 1e-12)
  expect_equal(a$p_h1, b$p_h1, tolerance = 1e-12)
  expect_equal(a$p_h2, b$p_h2, tolerance = 1e-12)
})

test_that("exact compose vs lr_dist_compose_R on synthetic features", {
  d1 <- list(log10_lr = c(-1.3, 0.4, 2.1),
             p_h1 = c(0.2, 0.5, 0.3), p_h2 = c(0.4, 0.4, 0.2))
  d2 <- list(log10_lr = c(0.1, 0.9),
             p_h1 = c(0.55, 0.45), p_h2 = c(0.25, 0.75))
  d3 <- list(log10_lr = c(-0.2, 0.6, 1.0, 1.7),
             p_h1 = c(0.1, 0.2, 0.3, 0.4), p_h2 = c(0.4, 0.3, 0.2, 0.1))
  ref <- mispitools:::lr_dist_compose_R(list(as_df(d1), as_df(d2), as_df(d3)))
  cpp <- mispitools:::cpp_lr_dist_compose(list(d1, d2, d3))
  expect_dist_equal(cpp, ref, info = "3 synthetic features")
})

test_that("method must be 'exact' or 'grid'", {
  d <- list(log10_lr = 0, p_h1 = 1, p_h2 = 1)
  expect_error(
    mispitools:::cpp_lr_dist_compose(list(d), method = "fft"),
    "'exact' or 'grid'"
  )
})

# ---------------------------------------------------------------------------
# Grid heuristic: mass + mean preserved exactly; shape converges to exact.
# ---------------------------------------------------------------------------

test_that("grid compose preserves total mass and the composed mean", {
  d1 <- list(log10_lr = c(-1.3, 0.4, 2.1),
             p_h1 = c(0.2, 0.5, 0.3), p_h2 = c(0.4, 0.4, 0.2))
  d2 <- list(log10_lr = c(0.1, 0.9),
             p_h1 = c(0.55, 0.45), p_h2 = c(0.25, 0.75))
  ex <- mispitools:::cpp_lr_dist_compose(list(d1, d2), method = "exact")
  gr <- mispitools:::cpp_lr_dist_compose(list(d1, d2), method = "grid",
                                         grid_points = 256)
  expect_equal(sum(gr$p_h1), 1, tolerance = 1e-12)
  expect_equal(sum(gr$p_h2), 1, tolerance = 1e-12)
  # Linear splitting preserves the first moment exactly.
  mean_ex <- sum(ex$p_h1 * ex$log10_lr)
  mean_gr <- sum(gr$p_h1 * gr$log10_lr)
  expect_equal(mean_gr, mean_ex, tolerance = 1e-9)
})

test_that("grid converges to exact as grid_points grows (1-Wasserstein gap)", {
  d1 <- list(log10_lr = c(-1.3, 0.4, 2.1),
             p_h1 = c(0.2, 0.5, 0.3), p_h2 = c(0.4, 0.4, 0.2))
  d2 <- list(log10_lr = c(0.1, 0.9),
             p_h1 = c(0.55, 0.45), p_h2 = c(0.25, 0.75))
  ex <- mispitools:::cpp_lr_dist_compose(list(d1, d2), method = "exact")
  # 1-Wasserstein distance = ∫ |CDF_grid - CDF_exact|. Both CDFs are
  # right-continuous step functions; evaluate on a fine common axis.
  cdf_at <- function(d, xs) {
    o <- order(d$log10_lr)
    stepfun(d$log10_lr[o], c(0, cumsum(d$p_h1[o])))(xs)
  }
  w1 <- function(g) {
    xs <- seq(min(ex$log10_lr, g$log10_lr),
              max(ex$log10_lr, g$log10_lr), length.out = 4000)
    sum(abs(cdf_at(g, xs) - cdf_at(ex, xs))) * (xs[2] - xs[1])
  }
  gr_coarse <- mispitools:::cpp_lr_dist_compose(list(d1, d2),
                                                method = "grid",
                                                grid_points = 16)
  gr_fine <- mispitools:::cpp_lr_dist_compose(list(d1, d2),
                                              method = "grid",
                                              grid_points = 1024)
  expect_lt(w1(gr_fine), w1(gr_coarse))
})

test_that("grid method rejects infinite supports", {
  d_inf <- list(log10_lr = c(-Inf, 0.5), p_h1 = c(0, 1), p_h2 = c(0.3, 0.7))
  expect_error(
    mispitools:::cpp_lr_dist_compose(list(d_inf), method = "grid"),
    "finite supports"
  )
})

# ---------------------------------------------------------------------------
# Cross-check vs R-ref on real per-marker joints + self-consistency.
# ---------------------------------------------------------------------------

test_that("exact compose matches R-ref over 3 Argentina markers (mut=equal)", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  mk <- function(fr) {
    marker_model(ped, "M", fr,
                 mutation = list(model = "equal", rate = 0.005))
  }
  m1 <- mk(c(a = 0.4, b = 0.6))
  m2 <- mk(c(a = 0.3, b = 0.5, c = 0.2))
  m3 <- mk(c(a = 0.1, b = 0.2, c = 0.3, d = 0.4))
  dl <- list(lrdist_of(m1), lrdist_of(m2), lrdist_of(m3))
  ref <- mispitools:::lr_dist_compose_R(lapply(dl, as_df))
  cpp <- mispitools:::cpp_lr_dist_compose(dl)
  expect_dist_equal(cpp, ref, info = "3 markers mut=equal")
  expect_equal(sum(cpp$p_h1), 1, tolerance = 1e-12)
  expect_equal(sum(cpp$p_h2), 1, tolerance = 1e-12)
})

test_that("composed mean = sum of per-marker E[log10 LR] (SCOUT_DNAtools)", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  mk <- function(fr) {
    marker_model(ped, "M", fr,
                 mutation = list(model = "equal", rate = 0.005))
  }
  models <- list(mk(c(a = 0.4, b = 0.6)),
                 mk(c(a = 0.3, b = 0.5, c = 0.2)),
                 mk(c(a = 0.1, b = 0.2, c = 0.3, d = 0.4)))
  dl <- lapply(models, lrdist_of)
  cpp <- mispitools:::cpp_lr_dist_compose(dl)
  kls <- vapply(models,
                function(m) mispitools:::per_marker_kl_R(m)$e_log10_lr_h1,
                numeric(1))
  expect_equal(sum(cpp$p_h1 * cpp$log10_lr), sum(kls), tolerance = 1e-10)
  kls2 <- vapply(models,
                 function(m) mispitools:::per_marker_kl_R(m)$e_log10_lr_h2,
                 numeric(1))
  expect_equal(sum(cpp$p_h2 * cpp$log10_lr), sum(kls2), tolerance = 1e-10)
})

test_that("exact compose propagates -Inf atoms (mut=none) to the total", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(2)
  m1 <- marker_model(ped, "M", c(a = 0.4, b = 0.6))            # mut=none
  m2 <- marker_model(ped, "M", c(a = 0.3, b = 0.5, c = 0.2),
                     mutation = list(model = "equal", rate = 0.005))
  dl <- list(lrdist_of(m1), lrdist_of(m2))
  ref <- mispitools:::lr_dist_compose_R(lapply(dl, as_df))
  cpp <- mispitools:::cpp_lr_dist_compose(dl)
  expect_dist_equal(cpp, ref, info = "mut=none ⊕ mut=equal")
  expect_true(cpp$has_neg_inf)
  expect_identical(cpp$log10_lr[1], -Inf)
})
