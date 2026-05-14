## F1.4 — cpt_marker_joint_R() with mutation = "equal".
## Cross-check against an analytic trio formula:
##   With K alleles and rate R, the mutation matrix is M[i,i] = 1 - R,
##   M[i,j] = R/(K-1) for i != j. The per-parent transmission probability
##   of allele x from genotype g is T_mut[g, x] = sum_y T[g, y] * M[y, x],
##   where T is the unmutated meiotic transmission.
##
## For unordered child genotype (c1, c2):
##   P(child | father f, mother m) = T_mut[f, c1] * T_mut[m, c2]
##                                   + (c1 != c2) * T_mut[f, c2] * T_mut[m, c1].

skip_if_no_pedtools <- function() {
  testthat::skip_if_not_installed("pedtools")
}

## Closed-form transmission probability from genotype g (a/b in unordered
## form) to allele x, under equal-rate mutation with rate R and K alleles.
## When K = 1, returns NA — undefined, but caller has K >= 2 by validation.
T_mut_equal <- function(g, x, R, K) {
  if (g[1L] == g[2L]) {
    if (g[1L] == x) 1 - R else R / (K - 1L)
  } else {
    half_left  <- if (g[1L] == x) 1 - R else R / (K - 1L)
    half_right <- if (g[2L] == x) 1 - R else R / (K - 1L)
    0.5 * half_left + 0.5 * half_right
  }
}

## Closed-form P(child genotype c | father f, mother m) under equal-rate
## mutation. c, f, m are character vectors of length 2 (alleles).
analytic_child_cond <- function(c, f, m, R, K) {
  pf1 <- T_mut_equal(f, c[1L], R, K)
  pm2 <- T_mut_equal(m, c[2L], R, K)
  if (c[1L] == c[2L]) {
    pf1 * pm2
  } else {
    pf2 <- T_mut_equal(f, c[2L], R, K)
    pm1 <- T_mut_equal(m, c[1L], R, K)
    pf1 * pm2 + pf2 * pm1
  }
}

## Parse "x/y" -> c("x","y").
parse_geno <- function(g) strsplit(g, "/", fixed = TRUE)[[1L]]

## HWE prior for an unordered genotype label "x/y" given named freq vector.
hwe_prob <- function(g, freqs) {
  ab <- parse_geno(g)
  pa <- freqs[[ab[1L]]]; pb <- freqs[[ab[2L]]]
  if (ab[1L] == ab[2L]) pa^2 else 2 * pa * pb
}

trio_model_eq <- function(freqs, R) {
  marker_model(pedtools::nuclearPed(1), "M1", freqs,
               mutation = list(model = "equal", rate = R))
}

# ---------------------------------------------------------------------------
# 1) Sanity / structural checks under equal-rate mutation
# ---------------------------------------------------------------------------

test_that("P_H1 and P_H2 each sum to 1 under equal-rate mutation (K=2)", {
  skip_if_no_pedtools()
  out <- mispitools:::cpt_marker_joint_R(
    trio_model_eq(c("a" = 0.4, "b" = 0.6), R = 0.01))
  expect_equal(sum(out$P_H1), 1, tolerance = 1e-12)
  expect_equal(sum(out$P_H2), 1, tolerance = 1e-12)
  expect_true(all(out$P_H1 >= 0))
  expect_true(all(out$P_H2 >= 0))
})

test_that("P_H1 and P_H2 each sum to 1 under equal-rate mutation (K=4)", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.25, "b" = 0.25, "c" = 0.25, "d" = 0.25)
  out <- mispitools:::cpt_marker_joint_R(trio_model_eq(freqs, R = 0.005))
  expect_equal(sum(out$P_H1), 1, tolerance = 1e-12)
  expect_equal(sum(out$P_H2), 1, tolerance = 1e-12)
})

test_that("equal-rate mutation: support is full (no Mendelian-zero rows)", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.6)
  R <- 0.02
  out <- mispitools:::cpt_marker_joint_R(trio_model_eq(freqs, R))
  ## With R > 0 every parent x parent x child combination has positive
  ## P_H1 mass — mutation lifts the Mendelian zeros that exist under
  ## mutation = "none".
  expect_true(all(out$P_H1 > 0))
})

# ---------------------------------------------------------------------------
# 2) Closed-form trio: AA x AA -> child, K=2
# ---------------------------------------------------------------------------

test_that("analytic trio: parents (a/a, a/a), child distribution, K=2", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.6)
  R <- 0.01
  out <- mispitools:::cpt_marker_joint_R(trio_model_eq(freqs, R))

  ## Father transmits a w.p. 1-R, b w.p. R. Mother same.
  ## P(child=a/a) = (1-R)^2, P(child=a/b) = 2*R*(1-R), P(child=b/b) = R^2.
  pP <- freqs[["a"]]^2 * freqs[["a"]]^2

  r_aa <- out[out$`1` == "a/a" & out$`2` == "a/a" & out$`3` == "a/a", ]
  r_ab <- out[out$`1` == "a/a" & out$`2` == "a/a" & out$`3` == "a/b", ]
  r_bb <- out[out$`1` == "a/a" & out$`2` == "a/a" & out$`3` == "b/b", ]

  expect_equal(r_aa$P_H1, pP * (1 - R)^2,        tolerance = 1e-12)
  expect_equal(r_ab$P_H1, pP * 2 * R * (1 - R),  tolerance = 1e-12)
  expect_equal(r_bb$P_H1, pP * R^2,              tolerance = 1e-12)
  expect_equal(r_aa$P_H1 + r_ab$P_H1 + r_bb$P_H1, pP, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 3) Closed-form trio: AA x BB -> child, K=2
# ---------------------------------------------------------------------------

test_that("analytic trio: parents (a/a, b/b), child distribution, K=2", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.3, "b" = 0.7)
  R <- 0.05
  out <- mispitools:::cpt_marker_joint_R(trio_model_eq(freqs, R))

  ## Father (a/a) -> a w.p. 1-R, b w.p. R.
  ## Mother (b/b) -> b w.p. 1-R, a w.p. R.
  pP <- freqs[["a"]]^2 * freqs[["b"]]^2
  p_aa <- (1 - R) * R
  p_ab <- (1 - R)^2 + R^2
  p_bb <- R * (1 - R)

  r_aa <- out[out$`1` == "a/a" & out$`2` == "b/b" & out$`3` == "a/a", ]
  r_ab <- out[out$`1` == "a/a" & out$`2` == "b/b" & out$`3` == "a/b", ]
  r_bb <- out[out$`1` == "a/a" & out$`2` == "b/b" & out$`3` == "b/b", ]

  expect_equal(r_aa$P_H1, pP * p_aa, tolerance = 1e-12)
  expect_equal(r_ab$P_H1, pP * p_ab, tolerance = 1e-12)
  expect_equal(r_bb$P_H1, pP * p_bb, tolerance = 1e-12)
  expect_equal(p_aa + p_ab + p_bb, 1, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 4) Closed-form trio: AB x AB -> child, K=3
# ---------------------------------------------------------------------------

test_that("analytic trio: parents (a/b, a/b), child distribution, K=3", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.3, "b" = 0.5, "c" = 0.2)
  R <- 0.02
  out <- mispitools:::cpt_marker_joint_R(trio_model_eq(freqs, R))

  pP <- (2 * freqs[["a"]] * freqs[["b"]])^2
  parents <- list(f = c("a", "b"), m = c("a", "b"))
  K <- length(freqs)
  for (cg in c("a/a", "a/b", "a/c", "b/b", "b/c", "c/c")) {
    cv <- parse_geno(cg)
    expected_cond <- analytic_child_cond(cv, parents$f, parents$m, R, K)
    row <- out[out$`1` == "a/b" & out$`2` == "a/b" & out$`3` == cg, ]
    expect_equal(row$P_H1, pP * expected_cond, tolerance = 1e-12,
                 info = paste("child =", cg))
  }
})

# ---------------------------------------------------------------------------
# 5) Full closed-form check over every (father, mother, child) triple, K=3
# ---------------------------------------------------------------------------

test_that("analytic trio: full P_H1 matches closed form over all triples, K=3", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.3, "b" = 0.5, "c" = 0.2)
  R <- 0.015
  K <- length(freqs)
  out <- mispitools:::cpt_marker_joint_R(trio_model_eq(freqs, R))

  genos <- c("a/a", "a/b", "a/c", "b/b", "b/c", "c/c")
  for (fg in genos) for (mg in genos) for (cg in genos) {
    fv <- parse_geno(fg); mv <- parse_geno(mg); cv <- parse_geno(cg)
    expected <- hwe_prob(fg, freqs) * hwe_prob(mg, freqs) *
                analytic_child_cond(cv, fv, mv, R, K)
    row <- out[out$`1` == fg & out$`2` == mg & out$`3` == cg, ]
    expect_equal(row$P_H1, expected, tolerance = 1e-12,
                 info = paste(fg, "x", mg, "->", cg))
  }
})

# ---------------------------------------------------------------------------
# 6) P_H2 under mutation: factorization P_H2 = marg_parents * HWE(child)
# ---------------------------------------------------------------------------

test_that("P_H2 = sum_parents(P_H1) * HWE(child) under equal-rate mutation", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.3, "b" = 0.5, "c" = 0.2)
  R <- 0.01
  out <- mispitools:::cpt_marker_joint_R(trio_model_eq(freqs, R))

  ## POI = child (default for nuclearPed(1)). H2 = marg_parents(H1) * HWE(g_POI).
  parent_pairs <- unique(out[, c("1", "2")])
  for (i in seq_len(nrow(parent_pairs))) {
    rows <- out$`1` == parent_pairs$`1`[i] & out$`2` == parent_pairs$`2`[i]
    marg <- sum(out$P_H1[rows])
    sub <- out[rows, , drop = FALSE]
    hwe_child <- vapply(sub$`3`, hwe_prob, numeric(1), freqs = freqs)
    expect_equal(unname(sub$P_H2), unname(marg * hwe_child),
                 tolerance = 1e-12)
  }
})

# ---------------------------------------------------------------------------
# 7) Continuity: R -> 0 reproduces mutation = "none" on shared support
# ---------------------------------------------------------------------------

test_that("equal-rate mutation R -> 0 matches mutation=none on shared support", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.3, "b" = 0.5, "c" = 0.2)
  mm_none <- marker_model(pedtools::nuclearPed(1), "M1", freqs)
  mm_tiny <- trio_model_eq(freqs, R = 1e-12)
  out_none <- mispitools:::cpt_marker_joint_R(mm_none)
  out_tiny <- mispitools:::cpt_marker_joint_R(mm_tiny)

  key_none <- do.call(paste, c(out_none[, c("1", "2", "3")], list(sep = "_")))
  key_tiny <- do.call(paste, c(out_tiny[, c("1", "2", "3")], list(sep = "_")))
  for (k in key_none) {
    pn <- out_none$P_H1[match(k, key_none)]
    pt <- out_tiny$P_H1[match(k, key_tiny)]
    expect_equal(pn, pt, tolerance = 1e-10, info = paste("row", k))
  }
})

# ---------------------------------------------------------------------------
# 8) Half-sib pedigree + equal-rate mutation: structural validity
# ---------------------------------------------------------------------------

test_that("half-sib pedigree with equal-rate mutation: joint sums to 1", {
  skip_if_no_pedtools()
  ped <- pedtools::halfSibPed()
  freqs <- c("a" = 0.3, "b" = 0.4, "c" = 0.3)
  mm <- marker_model(ped, "M1", freqs,
                     mutation = list(model = "equal", rate = 0.01))
  out <- mispitools:::cpt_marker_joint_R(mm)
  expect_equal(sum(out$P_H1), 1, tolerance = 1e-10)
  expect_equal(sum(out$P_H2), 1, tolerance = 1e-10)
  expect_true(all(out$P_H1 > 0))
  expect_true(all(out$P_H2 > 0))
})

# ---------------------------------------------------------------------------
# 9) Marginal allele frequency under H1 is preserved by equal-rate mutation
# ---------------------------------------------------------------------------

test_that("marginal allele frequency of the child equals input freqs", {
  skip_if_no_pedtools()
  ## Equal-rate mutation is symmetric: a population at HWE on allele
  ## frequencies p stays at marginal p after one generation of mutated
  ## transmission, because sum_a p_a * M[a,x] = p_x*(1-R) + (1-p_x)*R/(K-1)
  ## equals p_x only when M satisfies stationary detailed balance with p.
  ## For uniform p, this holds exactly. Use uniform freqs to test.
  freqs <- c("a" = 1/3, "b" = 1/3, "c" = 1/3)
  R <- 0.03
  out <- mispitools:::cpt_marker_joint_R(trio_model_eq(freqs, R))

  ## Marginal P(child genotype) = sum over parents of P_H1.
  child_marg <- as.numeric(tapply(out$P_H1, out$`3`, sum))
  ## HWE on (uniform) freqs:
  expected <- as.numeric(vapply(sort(unique(out$`3`)), hwe_prob,
                                numeric(1), freqs = freqs))
  expect_equal(child_marg, expected, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# 10) Mutation matrix builder is well-formed
# ---------------------------------------------------------------------------

test_that("mutation_matrix_R(equal) is stochastic with correct entries", {
  for (K in c(2L, 3L, 5L, 10L)) for (R in c(0, 1e-6, 1e-3, 0.1, 0.5)) {
    M <- mispitools:::mutation_matrix_R(list(model = "equal", rate = R), K)
    expect_equal(dim(M), c(K, K))
    expect_equal(rowSums(M), rep(1, K), tolerance = 1e-12)
    expect_equal(colSums(M), rep(1, K), tolerance = 1e-12)
    expect_equal(diag(M), rep(1 - R, K),  tolerance = 1e-12)
    off <- M[lower.tri(M) | upper.tri(M)]
    expect_equal(unique(off), R / (K - 1), tolerance = 1e-12)
  }
})

test_that("mutation_matrix_R(none) is the identity", {
  for (K in c(2L, 3L, 7L)) {
    M <- mispitools:::mutation_matrix_R(list(model = "none", rate = 0), K)
    expect_equal(M, diag(K))
  }
})
