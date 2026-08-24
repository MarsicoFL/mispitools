## F4.4 — lr_distribution() public API + quantile.lr_dist() method.
##
## Oracle chain: lr_distribution() composes per-marker LR distributions
## through cpp_lr_dist_compose; the R-reference is lr_dist_compose_R()
## over per_marker_lr_dist_R() outputs (validated in F4.1/F4.2). Cross-
## checks pin the public layer bit-for-bit to that reference, and the
## composed mean to the F3 KL engine (Sigma e_log10_lr).

skip_if_not_installed("pedtools")

mm <- function(ped, id, freqs, mut = list(model = "none", rate = 0)) {
  marker_model(ped, id, freqs, mutation = mut)
}

ref_compose <- function(models, poi = NULL) {
  ds <- lapply(models, function(m)
    mispitools:::per_marker_lr_dist_R(m, poi = poi, aggregate = TRUE))
  mispitools:::lr_dist_compose_R(ds)
}

# ---------------------------------------------------------------------------
# Structure & attributes
# ---------------------------------------------------------------------------

test_that("lr_distribution returns an lr_dist with profile attributes", {
  ped <- ped_parent_child()
  models <- list(
    A = mm(ped, "A", c("1" = 0.4, "2" = 0.6),
           list(model = "equal", rate = 0.005)),
    B = mm(ped, "B", c("1" = 0.2, "2" = 0.3, "3" = 0.5),
           list(model = "equal", rate = 0.005))
  )
  d <- lr_distribution(models)
  expect_s3_class(d, "lr_dist")
  expect_s3_class(d, "data.frame")
  expect_named(d, c("log10_lr", "p_h1", "p_h2"))
  expect_equal(attr(d, "markers"), c("A", "B"))
  expect_equal(attr(d, "n_markers"), 2L)
  expect_equal(attr(d, "method"), "exact")
  expect_false(is.null(attr(d, "poi")))
  expect_equal(sum(d$p_h1), 1, tolerance = 1e-12)
  expect_equal(sum(d$p_h2), 1, tolerance = 1e-12)
})

test_that("a single marker_model is accepted and equals a one-element list", {
  ped <- ped_parent_child()
  m <- mm(ped, "M", c("1" = 0.3, "2" = 0.7),
          list(model = "equal", rate = 0.01))
  d_single <- lr_distribution(m)
  d_list <- lr_distribution(list(m))
  expect_equal(as.data.frame(d_single), as.data.frame(d_list),
               tolerance = 1e-15)
  expect_equal(attr(d_single, "n_markers"), 1L)
})

# ---------------------------------------------------------------------------
# Cross-check vs the R reference composition
# ---------------------------------------------------------------------------

test_that("single-model composition equals per_marker_lr_dist_R", {
  ped <- ped_parent_child()
  m <- mm(ped, "M", c("1" = 0.3, "2" = 0.7),
          list(model = "equal", rate = 0.005))
  d <- lr_distribution(list(m))
  ref <- mispitools:::per_marker_lr_dist_R(m, aggregate = TRUE)
  expect_equal(d$log10_lr, ref$log10_lr, tolerance = 1e-12)
  expect_equal(d$p_h1, ref$p_h1, tolerance = 1e-12)
  expect_equal(d$p_h2, ref$p_h2, tolerance = 1e-12)
})

test_that("multi-marker exact composition matches lr_dist_compose_R", {
  peds <- list(ped_parent_child(), ped_half_sibs(),
               ped_grandparent_grandchild())
  for (ped in peds) {
    models <- list(
      mm(ped, "M1", c("1" = 0.4, "2" = 0.6),
         list(model = "equal", rate = 0.005)),
      mm(ped, "M2", c("1" = 0.2, "2" = 0.3, "3" = 0.5),
         list(model = "equal", rate = 0.01)),
      mm(ped, "M3", c("10" = 0.5, "11" = 0.5),
         list(model = "stepwise", rate = 0.004, ratio = 0.1))
    )
    d <- lr_distribution(models)
    ref <- ref_compose(models)
    expect_equal(nrow(d), nrow(ref))
    expect_equal(d$log10_lr, ref$log10_lr, tolerance = 1e-12)
    expect_equal(d$p_h1, ref$p_h1, tolerance = 1e-12)
    expect_equal(d$p_h2, ref$p_h2, tolerance = 1e-12)
  }
})

test_that("Argentina markers compose against the reference", {
  data("Argentina", package = "mispitools", envir = environment())
  ped <- ped_parent_child()
  fA <- top_k_freqs(Argentina, "D3S1358", 3L)
  fB <- top_k_freqs(Argentina, "D5S818", 4L)
  models <- list(
    D3S1358 = mm(ped, "D3S1358", fA, list(model = "equal", rate = 0.002)),
    D5S818  = mm(ped, "D5S818",  fB, list(model = "equal", rate = 0.002))
  )
  d <- lr_distribution(models)
  ref <- ref_compose(models)
  expect_equal(d$log10_lr, ref$log10_lr, tolerance = 1e-10)
  expect_equal(d$p_h1, ref$p_h1, tolerance = 1e-10)
  expect_equal(d$p_h2, ref$p_h2, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Self-consistency with the F3 KL engine
# ---------------------------------------------------------------------------

test_that("composed E[log10 LR | Hk] == sum of per-marker expectations", {
  ped <- ped_parent_child()
  models <- list(
    mm(ped, "M1", c("1" = 0.4, "2" = 0.6),
       list(model = "equal", rate = 0.005)),
    mm(ped, "M2", c("1" = 0.2, "2" = 0.3, "3" = 0.5),
       list(model = "equal", rate = 0.005)),
    mm(ped, "M3", c("1" = 0.5, "2" = 0.5),
       list(model = "equal", rate = 0.01))
  )
  d <- lr_distribution(models)
  s <- summary(d)
  kl <- per_marker_kl_profile(models)
  expect_equal(s$mean_h1, sum(kl$e_log10_lr_h1), tolerance = 1e-10)
  expect_equal(s$mean_h2, sum(kl$e_log10_lr_h2), tolerance = 1e-10)
})

test_that("summary() and plot() apply to the composed object", {
  ped <- ped_parent_child()
  models <- list(
    mm(ped, "M1", c("1" = 0.4, "2" = 0.6),
       list(model = "equal", rate = 0.005)),
    mm(ped, "M2", c("1" = 0.3, "2" = 0.7),
       list(model = "equal", rate = 0.005))
  )
  d <- lr_distribution(models)
  s <- summary(d)
  expect_s3_class(s, "summary.lr_dist")
  expect_true(is.finite(s$auc) && s$auc >= 0 && s$auc <= 1)
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(d), "ggplot")
})

# ---------------------------------------------------------------------------
# Infinite atoms (mutation = none) propagate through exact composition
# ---------------------------------------------------------------------------

test_that("mutation=none keeps -Inf atoms; grid rejects infinite support", {
  ped <- ped_parent_child()
  models <- list(
    mm(ped, "M1", c("1" = 0.5, "2" = 0.5)),
    mm(ped, "M2", c("1" = 0.4, "2" = 0.6))
  )
  d <- lr_distribution(models, method = "exact")
  expect_true(any(is.infinite(d$log10_lr)))
  expect_equal(sum(d$p_h1), 1, tolerance = 1e-12)
  expect_error(lr_distribution(models, method = "grid"))
})

test_that("grid method preserves mass and mean on finite support", {
  ped <- ped_parent_child()
  models <- list(
    mm(ped, "M1", c("1" = 0.4, "2" = 0.6),
       list(model = "equal", rate = 0.005)),
    mm(ped, "M2", c("1" = 0.2, "2" = 0.3, "3" = 0.5),
       list(model = "equal", rate = 0.005))
  )
  d_ex <- lr_distribution(models, method = "exact")
  d_gr <- lr_distribution(models, method = "grid", grid_points = 2048L)
  expect_equal(sum(d_gr$p_h1), 1, tolerance = 1e-9)
  mean_ex <- stats::weighted.mean(d_ex$log10_lr, d_ex$p_h1)
  mean_gr <- stats::weighted.mean(d_gr$log10_lr, d_gr$p_h1)
  expect_equal(mean_gr, mean_ex, tolerance = 1e-9)
  expect_equal(attr(d_gr, "method"), "grid")
})

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

test_that("lr_distribution validates its inputs", {
  ped <- ped_parent_child()
  m <- mm(ped, "M", c("1" = 0.5, "2" = 0.5),
          list(model = "equal", rate = 0.005))
  expect_error(lr_distribution(42), "marker_model")
  expect_error(lr_distribution(list()), "at least one")
  expect_error(lr_distribution(list(m, 1L)), "position")
  expect_error(lr_distribution(list(m), poi = c("a", "b")),
               "single non-empty")
  expect_error(lr_distribution(list(m), poi = ""), "single non-empty")
  expect_error(lr_distribution(list(m), grid_points = 1L),
               "grid_points")
  expect_error(lr_distribution(list(m), method = "fft"))
})

# ---------------------------------------------------------------------------
# quantile.lr_dist
# ---------------------------------------------------------------------------

test_that("quantile.lr_dist matches the reference, both hypotheses", {
  ped <- ped_parent_child()
  models <- list(
    mm(ped, "M1", c("1" = 0.4, "2" = 0.6),
       list(model = "equal", rate = 0.005)),
    mm(ped, "M2", c("1" = 0.2, "2" = 0.3, "3" = 0.5),
       list(model = "equal", rate = 0.005))
  )
  d <- lr_distribution(models)
  ref <- ref_compose(models)
  probs <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
  q1 <- quantile(d, probs = probs, under_h1 = TRUE)
  q2 <- quantile(d, probs = probs, under_h1 = FALSE)
  expect_equal(unname(q1),
               mispitools:::lr_dist_quantile_R(ref, probs, under_h1 = TRUE),
               tolerance = 1e-12)
  expect_equal(unname(q2),
               mispitools:::lr_dist_quantile_R(ref, probs, under_h1 = FALSE),
               tolerance = 1e-12)
  expect_equal(names(q1),
               c("0%", "10%", "25%", "50%", "75%", "90%", "100%"))
  expect_true(all(diff(q1) >= 0))
})

test_that("quantile.lr_dist default probs are the quartiles", {
  d <- as_lr_dist(data.frame(
    log10_lr = c(-1, 0, 2),
    p_h1 = c(0.1, 0.3, 0.6),
    p_h2 = c(0.6, 0.3, 0.1)
  ))
  q <- quantile(d)
  expect_length(q, 5L)
  expect_equal(names(q), c("0%", "25%", "50%", "75%", "100%"))
})

test_that("quantile.lr_dist validates probs and under_h1", {
  d <- as_lr_dist(data.frame(
    log10_lr = c(-1, 1), p_h1 = c(0.5, 0.5), p_h2 = c(0.5, 0.5)))
  expect_error(quantile(d, probs = c(-0.1, 0.5)), "\\[0, 1\\]")
  expect_error(quantile(d, probs = 1.2), "\\[0, 1\\]")
  expect_error(quantile(d, probs = "x"), "numeric")
  expect_error(quantile(d, under_h1 = NA), "single logical")
})

test_that("quantile.lr_dist errors with no support under the hypothesis", {
  d <- as_lr_dist(data.frame(
    log10_lr = c(-1, 1), p_h1 = c(0, 0), p_h2 = c(0.5, 0.5)))
  expect_error(quantile(d, under_h1 = TRUE))
  q <- quantile(d, probs = 0.5, under_h1 = FALSE)
  expect_equal(unname(q), -1)
})
