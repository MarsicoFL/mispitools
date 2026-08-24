## F1.5 — Stepwise mutation in the R reference engine.
##
## Two oracles for the per-locus mutation matrix:
##   * fbnet::getLocusCPT (the spec-mandated F1.5 oracle), reshaped to
##     a K x K matrix mij[Ap, node].
##   * pedmut::mutationMatrix(model = "stepwise", rate2 = 0), the
##     canonical Familias-style stepwise model on integer alleles.
##
## Joint CPT oracle: pedprobr::oneMarkerDistribution() with a stepwise
## mutation model attached to the marker via pedmut. Tested on the
## half-sib pedigree with D3S1358 trimmed to its top-3 alleles to keep
## the brute-force enumeration tractable (every off-diagonal entry of
## the stepwise matrix is positive, so the engine cannot prune).

skip_if_no_pedtools <- function() {
  testthat::skip_if_not_installed("pedtools")
}

extract_fbnet_mij <- function(freqs, marker_name, R, r) {
  bn <- list(
    alelFreq = stats::setNames(list(freqs), marker_name),
    mmodel = list(
      hombre = list(stepwise = c(R = R, R2 = 0, r = r)),
      mujer  = list(stepwise = c(R = R, R2 = 0, r = r))
    )
  )
  out <- fbnet::getLocusCPT(bn, marker_name)
  alleles <- names(freqs)
  K <- length(alleles)
  sub <- out[out$S == 1L & out$Am == alleles[1L], ]
  M <- matrix(0, nrow = K, ncol = K, dimnames = list(alleles, alleles))
  for (k in seq_len(nrow(sub))) {
    M[sub$Ap[k], sub$node[k]] <- sub$prob[k]
  }
  M
}

joint_max_abs_diff_stepwise <- function(ped, marker_name, freqs, R, r) {
  ids <- as.character(labels(ped))

  mut_spec <- list(model = "stepwise", rate = R, ratio = r)
  mm <- marker_model(ped, marker_name, freqs, mutation = mut_spec)
  tab <- mispitools:::cpt_marker_joint_R(mm)

  pedmut_mat <- pedmut::mutationMatrix(
    model = "stepwise", alleles = names(freqs),
    rate = R, rate2 = 0, range = r
  )
  m_obj <- pedtools::marker(ped, afreq = freqs, name = marker_name,
                            mutmod = pedmut_mat)
  ped_with <- pedtools::addMarkers(ped, m_obj)
  oracle <- pedprobr::oneMarkerDistribution(
    ped_with, ids = ids, marker = marker_name,
    output = "table", verbose = FALSE
  )

  key_m <- do.call(paste, c(tab[, ids, drop = FALSE], list(sep = "|")))
  key_o <- do.call(paste, c(oracle[, ids, drop = FALSE], list(sep = "|")))
  m_prob <- tapply(tab$P_H1, key_m, sum)
  o_prob <- stats::setNames(oracle$prob, key_o)

  all_keys <- union(names(m_prob), names(o_prob))
  m_vec <- unname(m_prob[all_keys]); m_vec[is.na(m_vec)] <- 0
  o_vec <- unname(o_prob[all_keys]); o_vec[is.na(o_vec)] <- 0
  max(abs(m_vec - o_vec))
}

test_that("stepwise matrix matches fbnet::getLocusCPT (D3S1358 top-4)", {
  skip_if_no_pedtools()
  testthat::skip_if_not_installed("fbnet")
  data(Argentina, package = "mispitools", envir = environment())
  freqs <- top_k_freqs(Argentina, "D3S1358", top_k = 4L)
  alleles <- names(freqs)
  R <- 0.005; r <- 0.1

  M_mine <- mispitools:::mutation_matrix_R(
    list(model = "stepwise", rate = R, ratio = r),
    K = length(alleles), alleles = alleles
  )
  dimnames(M_mine) <- list(alleles, alleles)

  M_fbnet <- extract_fbnet_mij(freqs, "D3S1358", R = R, r = r)

  expect_equal(unname(M_mine), unname(M_fbnet), tolerance = 1e-12)
  expect_equal(unname(rowSums(M_mine)), rep(1, length(alleles)),
               tolerance = 1e-14)
})

test_that("stepwise matrix matches pedmut canonical (rate2 = 0)", {
  testthat::skip_if_not_installed("pedmut")
  data(Argentina, package = "mispitools", envir = environment())
  freqs <- top_k_freqs(Argentina, "D3S1358", top_k = 4L)
  alleles <- names(freqs)
  R <- 0.005; r <- 0.1

  M_mine <- mispitools:::mutation_matrix_R(
    list(model = "stepwise", rate = R, ratio = r),
    K = length(alleles), alleles = alleles
  )
  M_pedmut <- pedmut::mutationMatrix(
    model = "stepwise", alleles = alleles,
    rate = R, rate2 = 0, range = r
  )
  M_pedmut_bare <- matrix(as.numeric(M_pedmut),
                          nrow = length(alleles), ncol = length(alleles))
  expect_equal(unname(M_mine), M_pedmut_bare, tolerance = 1e-12)
})

test_that("stepwise matrix: 3-allele AB×AB analytic check", {
  alleles <- c("10", "11", "12")
  R <- 0.01; r <- 0.5
  M <- mispitools:::mutation_matrix_R(
    list(model = "stepwise", rate = R, ratio = r),
    K = 3L, alleles = alleles
  )
  ## row 10 -> 11,12 weights r^1, r^2 = 0.5, 0.25; sum = 0.75
  expect_equal(M[1, 1], 1 - R, tolerance = 1e-15)
  expect_equal(M[1, 2], R * 0.5  / 0.75, tolerance = 1e-15)
  expect_equal(M[1, 3], R * 0.25 / 0.75, tolerance = 1e-15)
  ## row 11 -> 10,12 weights 0.5, 0.5; sum = 1.0
  expect_equal(M[2, 1], R * 0.5, tolerance = 1e-15)
  expect_equal(M[2, 3], R * 0.5, tolerance = 1e-15)
  expect_equal(M[2, 2], 1 - R, tolerance = 1e-15)
  ## row 12 mirrors row 10
  expect_equal(M[3, 3], 1 - R, tolerance = 1e-15)
  expect_equal(M[3, 2], R * 0.5  / 0.75, tolerance = 1e-15)
  expect_equal(M[3, 1], R * 0.25 / 0.75, tolerance = 1e-15)
  expect_equal(unname(rowSums(M)), c(1, 1, 1), tolerance = 1e-15)
})

test_that("stepwise joint CPT (half-sib + D3S1358) sums to 1 for H1 and H2", {
  skip_if_no_pedtools()
  data(Argentina, package = "mispitools", envir = environment())
  freqs <- top_k_freqs(Argentina, "D3S1358", top_k = 3L)

  ped <- ped_half_sibs()
  mm <- marker_model(ped, "D3S1358", freqs,
                     mutation = list(model = "stepwise",
                                     rate = 0.005, ratio = 0.1))
  tab <- mispitools:::cpt_marker_joint_R(mm)

  expect_equal(sum(tab$P_H1), 1, tolerance = 1e-12)
  expect_equal(sum(tab$P_H2), 1, tolerance = 1e-12)
  expect_true(all(tab$P_H1 >= 0))
  expect_true(all(tab$P_H2 >= 0))
})

test_that("stepwise joint CPT matches pedprobr+pedmut on half-sib + D3S1358", {
  skip_if_no_pedtools()
  testthat::skip_if_not_installed("pedmut")
  testthat::skip_if_not_installed("pedprobr")
  data(Argentina, package = "mispitools", envir = environment())
  freqs <- top_k_freqs(Argentina, "D3S1358", top_k = 3L)
  ped <- ped_half_sibs()

  d <- joint_max_abs_diff_stepwise(ped, "D3S1358", freqs,
                                   R = 0.005, r = 0.1)
  expect_lt(d, 1e-10)
})

test_that("stepwise -> no-mutation continuity as rate -> 0 (parent-child)", {
  skip_if_no_pedtools()
  data(Argentina, package = "mispitools", envir = environment())
  freqs <- top_k_freqs(Argentina, "D3S1358", top_k = 3L)
  ped <- pedtools::nuclearPed(1)

  mm_none <- marker_model(ped, "D3S1358", freqs)
  mm_step <- marker_model(ped, "D3S1358", freqs,
                          mutation = list(model = "stepwise",
                                          rate = 1e-9, ratio = 0.1))
  tab_none <- mispitools:::cpt_marker_joint_R(mm_none)
  tab_step <- mispitools:::cpt_marker_joint_R(mm_step)

  ids <- as.character(labels(ped))
  key_n <- do.call(paste, c(tab_none[, ids, drop = FALSE], list(sep = "|")))
  key_s <- do.call(paste, c(tab_step[, ids, drop = FALSE], list(sep = "|")))
  pn <- stats::setNames(tab_none$P_H1, key_n)
  ps <- stats::setNames(tab_step$P_H1, key_s)
  all_keys <- union(names(pn), names(ps))
  vn <- unname(pn[all_keys]); vn[is.na(vn)] <- 0
  vs <- unname(ps[all_keys]); vs[is.na(vs)] <- 0
  expect_lt(max(abs(vn - vs)), 1e-7)
})

test_that("stepwise without ratio errors with informative message", {
  expect_error(
    mispitools:::mutation_matrix_R(
      list(model = "stepwise", rate = 1e-3),
      K = 3L, alleles = c("10", "11", "12")
    ),
    "ratio"
  )
})

test_that("stepwise with non-numeric allele labels errors", {
  expect_error(
    mispitools:::mutation_matrix_R(
      list(model = "stepwise", rate = 1e-3, ratio = 0.1),
      K = 3L, alleles = c("A", "B", "C")
    ),
    "numeric allele labels"
  )
})

test_that("stepwise via cpt_marker_joint_R requires alleles plumbing", {
  skip_if_no_pedtools()
  ## Non-numeric allele labels propagate from the model into the engine
  ## and trigger the matrix-builder validation.
  freqs <- c("A" = 0.4, "B" = 0.6)
  ped <- pedtools::nuclearPed(1)
  mm <- marker_model(ped, "M1", freqs,
                     mutation = list(model = "stepwise",
                                     rate = 1e-3, ratio = 0.1))
  expect_error(mispitools:::cpt_marker_joint_R(mm),
               "numeric allele labels")
})
