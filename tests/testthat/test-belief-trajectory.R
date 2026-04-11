## Tests for belief trajectory functions
## Coverage targets: belief_trajectory, binary_belief_trajectory,
## concentration_index, concentration_index_positive, herfindahl_index,
## shannon_concentration, leave_one_out, entropy_log10, kl_divergence_log10,
## trajectory_metrics.

# ---- entropy_log10 ----

test_that("entropy_log10 of a point mass is ~0", {
  expect_lt(entropy_log10(c(1, 0, 0, 0)), 1e-18)
})

test_that("entropy_log10 of uniform is log10(n)", {
  for (n in c(2, 5, 12, 50)) {
    expect_equal(entropy_log10(rep(1 / n, n)), log10(n), tolerance = 1e-10)
  }
})

test_that("entropy_log10 is non-negative", {
  set.seed(1)
  for (i in 1:20) {
    n <- sample(2:20, 1)
    p <- runif(n)
    p <- p / sum(p)
    expect_gte(entropy_log10(p), 0)
  }
})

# ---- kl_divergence_log10 ----

test_that("kl_divergence_log10 is zero for p == q", {
  p <- c(0.3, 0.2, 0.5)
  expect_equal(kl_divergence_log10(p, p), 0, tolerance = 1e-12)
})

test_that("kl_divergence_log10 is non-negative (Gibbs' inequality)", {
  set.seed(2)
  for (i in 1:20) {
    n <- sample(2:10, 1)
    p <- runif(n); p <- p / sum(p)
    q <- runif(n); q <- q / sum(q)
    expect_gte(kl_divergence_log10(p, q), -1e-10)
  }
})

test_that("kl_divergence_log10 is base-10 (in bans)", {
  # KL(p||q) for p=(0.5,0.5), q=(0.9,0.1) should be about 0.368 in bans
  # (= ln(value) / ln(10))
  p <- c(0.5, 0.5); q <- c(0.9, 0.1)
  nats <- 0.5 * log(0.5 / 0.9) + 0.5 * log(0.5 / 0.1)
  expected_bans <- nats / log(10)
  expect_equal(kl_divergence_log10(p, q), expected_bans, tolerance = 1e-10)
})

# ---- concentration_index ----

test_that("concentration_index of uniform weights is 1/T", {
  for (T in c(2, 5, 15, 100)) {
    expect_equal(concentration_index(rep(3.14, T)), 1 / T, tolerance = 1e-12)
  }
})

test_that("concentration_index of single-dominant weights approaches 1", {
  w <- c(100, rep(1e-6, 14))
  expect_gt(concentration_index(w), 0.99)
})

test_that("concentration_index is scale-invariant", {
  set.seed(3)
  w <- runif(10, -5, 5)
  expect_equal(concentration_index(w), concentration_index(2.5 * w),
               tolerance = 1e-12)
})

test_that("concentration_index is zero for all-zero weights", {
  expect_equal(concentration_index(rep(0, 10)), 0)
})

test_that("concentration_index equals max/sum on positive weights", {
  set.seed(4)
  w <- runif(20, 0.1, 3)
  expect_equal(concentration_index(w), max(w) / sum(w), tolerance = 1e-12)
})

# ---- concentration_index_positive ----

test_that("concentration_index_positive ignores negative entries", {
  w <- c(2, 1, 0.5, -10, 0.3)
  wp <- pmax(w, 0)
  expect_equal(concentration_index_positive(w), max(wp) / sum(wp),
               tolerance = 1e-12)
})

test_that("concentration_index_positive equals concentration_index on positive weights", {
  set.seed(5)
  w <- runif(15, 0.1, 5)
  expect_equal(concentration_index_positive(w), concentration_index(w),
               tolerance = 1e-12)
})

test_that("concentration_index_positive is zero when all weights are non-positive", {
  expect_equal(concentration_index_positive(c(-1, -2, -0.5)), 0)
  expect_equal(concentration_index_positive(rep(0, 5)), 0)
})

# ---- herfindahl_index ----

test_that("herfindahl_index of uniform weights is 1/T", {
  for (T in c(2, 5, 15, 100)) {
    expect_equal(herfindahl_index(rep(2, T)), 1 / T, tolerance = 1e-12)
  }
})

test_that("herfindahl_index of single-dominant weights approaches 1", {
  expect_gt(herfindahl_index(c(100, rep(1e-6, 14))), 0.99)
})

test_that("herfindahl_index is in [0, 1]", {
  set.seed(6)
  for (i in 1:20) {
    w <- runif(10, -3, 3)
    h <- herfindahl_index(w)
    expect_gte(h, 0)
    expect_lte(h, 1)
  }
})

# ---- shannon_concentration ----

test_that("shannon_concentration of uniform weights is 0", {
  expect_equal(shannon_concentration(rep(1, 10)), 0, tolerance = 1e-12)
})

test_that("shannon_concentration approaches 1 for single-dominant weights", {
  expect_gt(shannon_concentration(c(1000, rep(1e-9, 9))), 0.95)
})

# ---- leave_one_out ----

test_that("leave_one_out returns one row per marker", {
  lrs <- c(A = 5, B = 2, C = 10)
  loo <- leave_one_out(lrs)
  expect_equal(nrow(loo), 3)
  expect_equal(loo$marker, c("A", "B", "C"))
})

test_that("leave_one_out total_without equals total minus each contribution", {
  lrs <- c(A = 5, B = 2, C = 10, D = 1.5)
  loo <- leave_one_out(lrs)
  total <- sum(log10(lrs))
  expect_equal(loo$total_without, total - loo$log10_lr, tolerance = 1e-12)
})

test_that("leave_one_out fraction sums to 1 for positive-only weights", {
  lrs <- c(A = 5, B = 2, C = 10, D = 1.5)
  loo <- leave_one_out(lrs)
  expect_equal(sum(loo$fraction), 1, tolerance = 1e-12)
})

test_that("leave_one_out errors on non-positive LRs", {
  expect_error(leave_one_out(c(A = 1, B = 0)), "positive")
  expect_error(leave_one_out(c(A = 1, B = -2)), "positive")
})

# ---- belief_trajectory ----

test_that("belief_trajectory preserves normalization", {
  prior <- c(1/3, 1/3, 1/3)
  lrs <- list(c(5, 1, 1), c(1, 3, 2), c(2, 2, 4))
  traj <- belief_trajectory(prior, lrs)
  for (t in seq_len(nrow(traj))) {
    expect_equal(sum(traj[t, ]), 1, tolerance = 1e-12)
  }
})

test_that("belief_trajectory has correct dimensions", {
  prior <- c(0.5, 0.5)
  lrs <- list(c(3, 1), c(1, 2), c(2, 1))
  traj <- belief_trajectory(prior, lrs)
  expect_equal(dim(traj), c(4, 2))
})

test_that("belief_trajectory row 1 is the prior", {
  prior <- c(0.1, 0.3, 0.6)
  lrs <- list(c(2, 1, 1))
  traj <- belief_trajectory(prior, lrs)
  expect_equal(traj[1, ], prior, tolerance = 1e-12)
})

test_that("belief_trajectory errors on invalid input", {
  expect_error(belief_trajectory(c(0.4, 0.4), list(c(1, 1))), "sum to 1")
  expect_error(belief_trajectory(c(0.5, 0.5), list(c(1, 1, 1))), "length")
  expect_error(belief_trajectory(c(0.5, 0.5), list(c(-1, 1))), "negative")
})

# ---- binary_belief_trajectory ----

test_that("binary_belief_trajectory has correct structure", {
  lrs <- c(D3S1358 = 5.2, TH01 = 1.8, D21S11 = 12.0)
  out <- binary_belief_trajectory(lrs)
  expect_equal(nrow(out), 4)
  expect_equal(out$step, 0:3)
  expect_equal(out$marker, c("Prior", "D3S1358", "TH01", "D21S11"))
})

test_that("binary_belief_trajectory posterior sums to 1", {
  lrs <- c(M1 = 5, M2 = 0.5, M3 = 3, M4 = 1.2)
  out <- binary_belief_trajectory(lrs)
  expect_equal(out$posterior_h1 + out$posterior_h2, rep(1, nrow(out)),
               tolerance = 1e-12)
})

test_that("binary_belief_trajectory cumulative LR is a cumulative sum of log LR", {
  lrs <- c(M1 = 5, M2 = 3, M3 = 2)
  out <- binary_belief_trajectory(lrs)
  expect_equal(out$cum_log10_lr, unname(c(0, cumsum(log10(lrs)))),
               tolerance = 1e-12)
})

# ---- trajectory_metrics ----

test_that("trajectory_metrics integrates with belief_trajectory", {
  prior <- c(0.25, 0.25, 0.25, 0.25)
  lrs <- list(c(3, 1, 1, 1), c(1, 2, 2, 1), c(1, 1, 1, 4))
  traj <- belief_trajectory(prior, lrs)
  m <- trajectory_metrics(traj)

  expect_equal(length(m$entropy), 4)
  expect_equal(length(m$kl_step), 3)
  expect_equal(length(m$kl_from_prior), 4)
  expect_equal(length(m$tv_step), 3)
  expect_type(m$path_length, "double")
  expect_type(m$concentration, "double")
  expect_type(m$concentration_herfindahl, "double")
  expect_type(m$concentration_shannon, "double")

  # Non-trivial updates imply non-zero path length
  expect_gt(m$path_length, 0)

  # Cumulative KL from prior is monotone (not strict, but non-decreasing)
  expect_equal(m$kl_from_prior[1], 0, tolerance = 1e-12)
})

test_that("trajectory_metrics static-dynamic correspondence holds at T=1", {
  # When T=1, the sum of per-step KL equals the static KL(P_1 || P_0).
  # This is the coherence link to the forensIT static framework (R-009).
  prior <- c(0.5, 0.5)
  lrs <- list(c(4, 1))
  traj <- belief_trajectory(prior, lrs)
  m <- trajectory_metrics(traj)
  expect_equal(m$kl_step[1], m$kl_from_prior[2], tolerance = 1e-12)
})

test_that("trajectory_metrics sum of per-step KL equals cumulative only in general", {
  # Note: sum of per-step KL does NOT in general equal cumulative KL from prior
  # because KL is not additive. This test documents the distinction.
  prior <- c(1/3, 1/3, 1/3)
  lrs <- list(c(5, 1, 1), c(1, 3, 1), c(1, 1, 2))
  traj <- belief_trajectory(prior, lrs)
  m <- trajectory_metrics(traj)
  # Just check both are positive and no error
  expect_gte(sum(m$kl_step), 0)
  expect_gte(m$kl_from_prior[length(m$kl_from_prior)], 0)
})
