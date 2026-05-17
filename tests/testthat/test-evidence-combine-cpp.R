## F6.5 — core::evidence_combine(): case-level combination of per-feature
## LR distributions. Two modes:
##   * independent  — exact convolution (delegates to lr_dist_compose);
##     must match cpp_lr_dist_compose bit-for-bit.
##   * markov_se     — Egeland-Marsico (2026) chain: H1 stays independent,
##     H2 follows per-adjacent-pair transition matrices.
## The legacy self-oracle (lr_sex / lr_age / ...) and the Egeland-Marsico
## reproduction script are the verifier's job in F6.6 / F6.7. Here we pin
## the kernel contract: exact arithmetic, the independence-equivalence
## invariant, mass preservation, and the error surface.

## Build a per-feature LrDist list (log10_lr = log10(p1) - log10(p2),
## matching per_marker_lr_dist's boundary convention).
mk_feature <- function(p1, p2) {
  lr <- log10(p1) - log10(p2)
  list(log10_lr = lr, p_h1 = p1, p_h2 = p2)
}

## Brute-force the markov_se joint for two binary features (exact oracle).
ms_two_feature <- function(fa, fb, Tm) {
  rows <- list()
  for (i in seq_along(fa$p_h1)) {
    for (j in seq_along(fb$p_h1)) {
      h1 <- fa$p_h1[i] * fb$p_h1[j]
      h2 <- fa$p_h2[i] * Tm[i, j]
      if (h1 <= 0 && h2 <= 0) next
      lr <- if (h1 > 0 && h2 > 0) log10(h1) - log10(h2)
            else if (h1 > 0) Inf else -Inf
      rows[[length(rows) + 1L]] <- c(lr, h1, h2)
    }
  }
  m <- do.call(rbind, rows)
  o <- order(m[, 1])
  m <- m[o, , drop = FALSE]
  key <- m[, 1]
  agg <- tapply(seq_len(nrow(m)), match(key, unique(key)), function(idx) {
    c(key[idx[1]], sum(m[idx, 2]), sum(m[idx, 3]))
  })
  a <- do.call(rbind, agg)
  list(log10_lr = unname(a[, 1]),
       p_h1     = unname(a[, 2]),
       p_h2     = unname(a[, 3]))
}

test_that("independent mode == cpp_lr_dist_compose bit-for-bit", {
  fa <- mk_feature(c(0.9, 0.1), c(0.5, 0.5))
  fb <- mk_feature(c(0.7, 0.2, 0.1), c(0.4, 0.4, 0.2))
  fc <- mk_feature(c(0.6, 0.4), c(0.3, 0.7))
  dists <- list(fa, fb, fc)

  ref <- mispitools:::cpp_lr_dist_compose(dists)
  got <- mispitools:::cpp_evidence_combine(dists, mode = "independent")

  expect_identical(got$log10_lr, ref$log10_lr)
  expect_identical(got$p_h1, ref$p_h1)
  expect_identical(got$p_h2, ref$p_h2)
  expect_identical(got$has_pos_inf, ref$has_pos_inf)
  expect_identical(got$has_neg_inf, ref$has_neg_inf)
})

test_that("markov_se with marginal transitions reduces to independent", {
  fa <- mk_feature(c(0.9, 0.1), c(0.5, 0.5))
  fb <- mk_feature(c(0.8, 0.2), c(0.4, 0.6))
  dists <- list(fa, fb)

  ## Independence chain: every row equals fb's H2 marginal.
  Tindep <- matrix(c(0.4, 0.6,
                     0.4, 0.6), nrow = 2, byrow = TRUE)

  indep <- mispitools:::cpp_evidence_combine(dists, mode = "independent")
  ms <- mispitools:::cpp_evidence_combine(
    dists, mode = "markov_se", transition = list(Tindep))

  expect_equal(ms$log10_lr, indep$log10_lr, tolerance = 1e-12)
  expect_equal(ms$p_h1, indep$p_h1, tolerance = 1e-12)
  expect_equal(ms$p_h2, indep$p_h2, tolerance = 1e-12)
  ## Mass preserved.
  expect_equal(sum(ms$p_h1), 1, tolerance = 1e-12)
  expect_equal(sum(ms$p_h2), 1, tolerance = 1e-12)
})

test_that("markov_se reproduces the exact two-feature joint", {
  fa <- mk_feature(c(0.9, 0.1), c(0.5, 0.5))
  fb <- mk_feature(c(0.8, 0.2), c(0.4, 0.6))
  Tm <- matrix(c(0.7, 0.3,
                 0.2, 0.8), nrow = 2, byrow = TRUE)

  exp <- ms_two_feature(fa, fb, Tm)
  got <- mispitools:::cpp_evidence_combine(
    list(fa, fb), mode = "markov_se", transition = list(Tm))

  expect_equal(got$log10_lr, exp$log10_lr, tolerance = 1e-12)
  expect_equal(got$p_h1, exp$p_h1, tolerance = 1e-12)
  expect_equal(got$p_h2, exp$p_h2, tolerance = 1e-12)
  expect_equal(sum(got$p_h1), 1, tolerance = 1e-12)
  expect_equal(sum(got$p_h2), 1, tolerance = 1e-12)

  ## A genuinely dependent chain must change the combined distribution
  ## relative to independence (the paper's central point).
  indep <- mispitools:::cpp_evidence_combine(list(fa, fb),
                                             mode = "independent")
  e_h1 <- function(d) sum(d$p_h1 * d$log10_lr)
  expect_false(isTRUE(all.equal(e_h1(got), e_h1(indep))))
})

test_that("single feature is a no-op pass-through", {
  fa <- mk_feature(c(0.7, 0.2, 0.1), c(0.3, 0.3, 0.4))
  got <- mispitools:::cpp_evidence_combine(
    list(fa), mode = "markov_se", transition = list())
  ref <- mispitools:::cpp_evidence_combine(list(fa), mode = "independent")
  expect_equal(got$log10_lr, ref$log10_lr, tolerance = 1e-12)
  expect_equal(got$p_h1, ref$p_h1, tolerance = 1e-12)
  expect_equal(got$p_h2, ref$p_h2, tolerance = 1e-12)
})

test_that("invalid inputs flow back as errors, never throw across border", {
  fa <- mk_feature(c(0.9, 0.1), c(0.5, 0.5))
  fb <- mk_feature(c(0.8, 0.2), c(0.4, 0.6))

  expect_error(
    mispitools:::cpp_evidence_combine(list(fa, fb), mode = "bogus"))

  ## Wrong number of transition matrices (need K-1 = 1).
  expect_error(
    mispitools:::cpp_evidence_combine(
      list(fa, fb), mode = "markov_se", transition = list()))

  ## Non-stochastic row.
  Tbad <- matrix(c(0.7, 0.7,
                   0.2, 0.8), nrow = 2, byrow = TRUE)
  expect_error(
    mispitools:::cpp_evidence_combine(
      list(fa, fb), mode = "markov_se", transition = list(Tbad)))

  ## Wrong transition dimensions.
  Tdim <- matrix(c(0.5, 0.3, 0.2), nrow = 1)
  expect_error(
    mispitools:::cpp_evidence_combine(
      list(fa, fb), mode = "markov_se", transition = list(Tdim)))

  ## Empty input under markov_se.
  expect_error(
    mispitools:::cpp_evidence_combine(list(), mode = "markov_se"))
})
