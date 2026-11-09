
skip_if_not_installed("pedtools")

mm_typed <- function(ped, typed = NULL, freqs = c(A = 0.3, B = 0.7),
                     mut = list(model = "none", rate = 0)) {
  marker_model(ped, "M", freqs, mutation = mut, typed = typed)
}

finite_mean_h1 <- function(d) {
  fin <- is.finite(d$log10_lr)
  sum(d$log10_lr[fin] * d$p_h1[fin])
}

test_that("typed accepts NULL and rejects malformed input", {
  ped <- ped_parent_child()
  expect_null(mm_typed(ped)$typed)
  expect_error(mm_typed(ped, typed = "99"), "absent from the pedigree")
  expect_error(mm_typed(ped, typed = character(0)), "at least one member")
  expect_error(mm_typed(ped, typed = 1), "character vector")
  expect_error(mm_typed(ped, typed = c("1", NA)), "character vector")
})

test_that("typed is stored in canonical pedigree order without duplicates", {
  ped <- ped_parent_child()
  expect_equal(mm_typed(ped, typed = c("3", "1"))$typed, c("1", "3"))
  expect_equal(mm_typed(ped, typed = c("1", "1", "3"))$typed, c("1", "3"))
})

test_that("print reports the typing configuration", {
  ped <- ped_parent_child()
  expect_output(print(mm_typed(ped)), "typed *: every member")
  expect_output(print(mm_typed(ped, typed = "1")), "typed *: 1 of 3")
})

test_that("typed = every member equals typed = NULL, atom for atom", {
  ped <- ped_parent_child()
  d_dense  <- lr_distribution(mm_typed(ped), poi = "3")
  d_peeled <- lr_distribution(mm_typed(ped, typed = c("1", "2", "3")),
                              poi = "3")
  expect_equal(nrow(d_peeled), nrow(d_dense))
  fin <- is.finite(d_dense$log10_lr)
  expect_equal(d_peeled$log10_lr[fin], d_dense$log10_lr[fin],
               tolerance = 1e-12)
  expect_equal(d_peeled$p_h1, d_dense$p_h1, tolerance = 1e-12)
  expect_equal(d_peeled$p_h2, d_dense$p_h2, tolerance = 1e-12)
})

test_that("the joint table reports only the members that carry data", {
  ped <- ped_grandparent_grandchild()
  j <- mispitools:::cpt_marker_joint_cpp_wrap(mm_typed(ped, typed = "2"),
                                              poi = "5")
  expect_named(j, c("2", "5", "P_H1", "P_H2"))
  expect_equal(sum(j$P_H1), 1, tolerance = 1e-12)
  expect_equal(sum(j$P_H2), 1, tolerance = 1e-12)
  expect_equal(attr(j, "typed"), "2")
})

ref_mean_h1 <- function(ped, freqs, refs, poi) {
  m <- pedtools::marker(ped, afreq = freqs)
  pedm <- pedtools::setMarkers(ped, m)
  h1 <- pedprobr::oneMarkerDistribution(pedm, ids = c(refs, poi),
                                        marker = 1, verbose = FALSE)
  others <- setdiff(as.character(labels(ped)), poi)
  sub <- subset(pedm, others)
  pr <- pedprobr::oneMarkerDistribution(sub, ids = refs, marker = 1,
                                        verbose = FALSE)
  sing <- pedtools::singleton(poi)
  sing <- pedtools::setMarkers(sing, pedtools::marker(sing, afreq = freqs))
  ppoi <- pedprobr::oneMarkerDistribution(sing, ids = poi, marker = 1,
                                          verbose = FALSE)
  h2 <- outer(as.numeric(pr), as.numeric(ppoi))
  dim(h2) <- c(dim(pr), length(ppoi))
  dimnames(h2) <- c(dimnames(pr), list(names(ppoi)))
  v1 <- as.numeric(h1)
  v2 <- as.numeric(do.call(`[`, c(list(h2), dimnames(h1))))
  ok <- v1 > 0 & v2 > 0
  sum(v1[ok] * (log10(v1[ok]) - log10(v2[ok])))
}

test_that("one typed reference on a five-member pedigree matches pedprobr", {
  skip_if_not_installed("pedprobr")
  ped <- ped_grandparent_grandchild()
  f <- c(A = 0.3, B = 0.7)
  d <- lr_distribution(mm_typed(ped, typed = "2", freqs = f), poi = "5")
  expect_equal(finite_mean_h1(d), ref_mean_h1(ped, f, "2", "5"),
               tolerance = 1e-10)
  expect_equal(sum(d$p_h1), 1, tolerance = 1e-12)
  expect_equal(sum(d$p_h2), 1, tolerance = 1e-12)
})

test_that("two typed references on the same pedigree match pedprobr", {
  skip_if_not_installed("pedprobr")
  ped <- ped_grandparent_grandchild()
  f <- c(A = 0.4, B = 0.6)
  d <- lr_distribution(mm_typed(ped, typed = c("1", "2"), freqs = f),
                       poi = "5")
  expect_equal(finite_mean_h1(d), ref_mean_h1(ped, f, c("1", "2"), "5"),
               tolerance = 1e-10)
})

test_that("more typed references never reduce the expected weight of evidence", {
  skip_if_not_installed("pedprobr")
  ped <- ped_grandparent_grandchild()
  f <- c(A = 0.3, B = 0.7)
  one <- finite_mean_h1(lr_distribution(mm_typed(ped, typed = "2", freqs = f),
                                        poi = "5"))
  two <- finite_mean_h1(lr_distribution(
    mm_typed(ped, typed = c("1", "2"), freqs = f), poi = "5"))
  expect_gt(two, one)
})

test_that("the KL batch and the per-model loop agree under typed", {
  ped <- ped_grandparent_grandchild()
  models <- list(
    A = mm_typed(ped, typed = "2", freqs = c(A = 0.3, B = 0.7),
                 mut = list(model = "equal", rate = 1e-3)),
    B = mm_typed(ped, typed = "2", freqs = c(A = 0.2, B = 0.3, C = 0.5),
                 mut = list(model = "equal", rate = 1e-3))
  )
  batch <- per_marker_kl_profile(models, poi = "5")
  loop <- do.call(rbind, lapply(models, per_marker_kl, poi = "5"))
  expect_equal(batch$e_log10_lr_h1, loop$e_log10_lr_h1, tolerance = 1e-12)
  expect_equal(batch$kl_h1_to_h2, loop$kl_h1_to_h2, tolerance = 1e-12)
})

test_that("models with different typing configurations leave the batch path", {
  ped <- ped_grandparent_grandchild()
  models <- list(
    A = mm_typed(ped, typed = "2", mut = list(model = "equal", rate = 1e-3)),
    B = mm_typed(ped, typed = c("1", "2"),
                 mut = list(model = "equal", rate = 1e-3))
  )
  expect_null(mispitools:::per_marker_kl_profile_batch(models, poi = "5"))
  out <- per_marker_kl_profile(models, poi = "5")
  expect_equal(nrow(out), 2L)
  expect_gt(out$e_log10_lr_h1[2], out$e_log10_lr_h1[1])
})

test_that("composition under typed stays additive in the expected weight", {
  ped <- ped_grandparent_grandchild()
  models <- list(
    A = mm_typed(ped, typed = "2", freqs = c(A = 0.3, B = 0.7),
                 mut = list(model = "equal", rate = 1e-3)),
    B = mm_typed(ped, typed = "2", freqs = c(A = 0.5, B = 0.5),
                 mut = list(model = "equal", rate = 1e-3))
  )
  d <- lr_distribution(models, poi = "5")
  kl <- per_marker_kl_profile(models, poi = "5")
  expect_equal(finite_mean_h1(d), sum(kl$e_log10_lr_h1), tolerance = 1e-10)
  expect_equal(sum(d$p_h1), 1, tolerance = 1e-12)
})
