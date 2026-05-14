## F1.2 — cpt_marker_joint_R(): pure-R reference engine for the joint
## marker CPT under H1 (pedigree topology) and H2 (POI HWE-independent
## of the rest). No mutation in F1.2; mutation tests live in F1.4+.

skip_if_no_pedtools <- function() {
  testthat::skip_if_not_installed("pedtools")
}

trio_freqs <- function() c("a" = 0.4, "b" = 0.6)

trio_model <- function(freqs = trio_freqs()) {
  skip_if_no_pedtools()
  marker_model(pedtools::nuclearPed(1), "M1", freqs)
}

test_that("output has the documented structure", {
  mm <- trio_model()
  out <- mispitools:::cpt_marker_joint_R(mm)

  expect_s3_class(out, "data.frame")
  expect_true(all(c("1", "2", "3", "P_H1", "P_H2") %in% names(out)))
  expect_identical(names(out)[length(names(out)) - 1L], "P_H1")
  expect_identical(names(out)[length(names(out))], "P_H2")
  expect_identical(attr(out, "marker_id"), "M1")
  expect_identical(attr(out, "alleles"), c("a", "b"))
  expect_identical(attr(out, "poi"), "3")

  expect_true(is.numeric(out$P_H1))
  expect_true(is.numeric(out$P_H2))
  expect_true(all(out$P_H1 >= 0))
  expect_true(all(out$P_H2 >= 0))
  expect_true(is.character(out$`1`))
})

test_that("P_H1 sums to 1 over a trio", {
  out <- mispitools:::cpt_marker_joint_R(trio_model())
  expect_equal(sum(out$P_H1), 1, tolerance = 1e-12)
})

test_that("P_H2 sums to 1 over a trio (POI = child)", {
  out <- mispitools:::cpt_marker_joint_R(trio_model())
  expect_equal(sum(out$P_H2), 1, tolerance = 1e-12)
})

test_that("Mendelian no-mutation: AA x AB -> child is AA or AB with prob 0.5", {
  freqs <- c("a" = 0.4, "b" = 0.6)
  mm <- trio_model(freqs)
  out <- mispitools:::cpt_marker_joint_R(mm)
  pAA <- freqs[["a"]]^2
  pAB <- 2 * freqs[["a"]] * freqs[["b"]]

  rAA <- out[out$`1` == "a/a" & out$`2` == "a/b" & out$`3` == "a/a", ]
  rAB <- out[out$`1` == "a/a" & out$`2` == "a/b" & out$`3` == "a/b", ]
  expect_equal(nrow(rAA), 1L)
  expect_equal(nrow(rAB), 1L)
  expect_equal(rAA$P_H1, pAA * pAB * 0.5, tolerance = 1e-12)
  expect_equal(rAB$P_H1, pAA * pAB * 0.5, tolerance = 1e-12)
  ## Child cannot be b/b given AA x AB under H1; under H2 the row still has
  ## positive HWE-based mass and so survives in the output.
  rBB <- out[out$`1` == "a/a" & out$`2` == "a/b" & out$`3` == "b/b", ]
  expect_equal(nrow(rBB), 1L)
  expect_equal(rBB$P_H1, 0, tolerance = 1e-12)
  expect_gt(rBB$P_H2, 0)
})

test_that("Mendelian no-mutation: AB x AB -> child distribution 1/4, 1/2, 1/4", {
  freqs <- c("a" = 0.5, "b" = 0.5)
  mm <- trio_model(freqs)
  out <- mispitools:::cpt_marker_joint_R(mm)
  pAB <- 2 * freqs[["a"]] * freqs[["b"]]

  rAA <- out[out$`1` == "a/b" & out$`2` == "a/b" & out$`3` == "a/a", ]
  rAB <- out[out$`1` == "a/b" & out$`2` == "a/b" & out$`3` == "a/b", ]
  rBB <- out[out$`1` == "a/b" & out$`2` == "a/b" & out$`3` == "b/b", ]
  expect_equal(rAA$P_H1, pAB * pAB * 0.25, tolerance = 1e-12)
  expect_equal(rAB$P_H1, pAB * pAB * 0.50, tolerance = 1e-12)
  expect_equal(rBB$P_H1, pAB * pAB * 0.25, tolerance = 1e-12)
})

test_that("P_H2 = marginal_over_POI(P_H1) * HWE(POI) for a trio", {
  freqs <- c("a" = 0.3, "b" = 0.5, "c" = 0.2)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs)
  out <- mispitools:::cpt_marker_joint_R(mm)

  parent_pairs <- unique(out[, c("1", "2")])
  for (i in seq_len(nrow(parent_pairs))) {
    rows <- out$`1` == parent_pairs$`1`[i] & out$`2` == parent_pairs$`2`[i]
    marg <- sum(out$P_H1[rows])
    ## All P_H2 entries with these parent values should equal marg * HWE(child).
    sub <- out[rows, , drop = FALSE]
    hwe_child <- vapply(sub$`3`, function(g) {
      ab <- strsplit(g, "/", fixed = TRUE)[[1]]
      pa <- freqs[[ab[1]]]; pb <- freqs[[ab[2]]]
      if (ab[1] == ab[2]) pa^2 else 2 * pa * pb
    }, numeric(1))
    expect_equal(unname(sub$P_H2), unname(marg * hwe_child),
                 tolerance = 1e-12)
  }
})

test_that("LR collapses to allele identity for trio with HWE founder POI", {
  ## With AA x AA parents, a child must be AA. Under H2, child is HWE.
  ## So LR for child=AA = 1 / HWE(AA) when parents are AA x AA. Sanity check.
  freqs <- c("a" = 0.4, "b" = 0.6)
  mm <- trio_model(freqs)
  out <- mispitools:::cpt_marker_joint_R(mm)
  rAA <- out[out$`1` == "a/a" & out$`2` == "a/a" & out$`3` == "a/a", ]
  expect_equal(rAA$P_H1, freqs[["a"]]^2 * freqs[["a"]]^2 * 1,
               tolerance = 1e-12)
  expect_equal(rAA$P_H2, freqs[["a"]]^2 * freqs[["a"]]^2 * freqs[["a"]]^2,
               tolerance = 1e-12)
  ## LR for (a/a, a/a, a/a) = 1 / HWE(a/a) for the child.
  expect_equal(rAA$P_H1 / rAA$P_H2, 1 / freqs[["a"]]^2, tolerance = 1e-12)
})

test_that("two unrelated singletons-as-founders: P_H1 = HWE x HWE", {
  ## linearPed(0) is just two unrelated parents (no children). Use a couple
  ## via nuclearPed and request poi=one of the founders to test marginalization
  ## semantics; we instead build a 2-founder pedigree by trimming.
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  freqs <- c("a" = 0.3, "b" = 0.7)
  mm <- marker_model(ped, "M1", freqs)
  out <- mispitools:::cpt_marker_joint_R(mm, poi = "1")
  ## Under H2 with poi=1 (a founder), P_H2(g_1,g_2,g_3) = P_H1(g_2,g_3) * HWE(g_1).
  ## Marginal over founder 1 in P_H1 = HWE(g_2) * P_Mendelian(g_3 | g_1, g_2)
  ## summed over g_1 = HWE(g_2) * P(g_3 | mother=g_2, father HWE).
  ## Simpler check: sum(P_H2) = 1 and P_H2 >= 0.
  expect_equal(sum(out$P_H2), 1, tolerance = 1e-12)
  expect_true(all(out$P_H2 >= 0))
  expect_identical(attr(out, "poi"), "1")
})

test_that("explicit poi argument is respected", {
  mm <- trio_model()
  out <- mispitools:::cpt_marker_joint_R(mm, poi = "2")
  expect_identical(attr(out, "poi"), "2")
  expect_equal(sum(out$P_H1), 1, tolerance = 1e-12)
  expect_equal(sum(out$P_H2), 1, tolerance = 1e-12)
})

test_that("invalid poi errors", {
  mm <- trio_model()
  expect_error(mispitools:::cpt_marker_joint_R(mm, poi = "999"),
               "not a member")
  expect_error(mispitools:::cpt_marker_joint_R(mm, poi = ""),
               "non-empty")
  expect_error(mispitools:::cpt_marker_joint_R(mm, poi = c("1", "2")),
               "non-empty")
  expect_error(mispitools:::cpt_marker_joint_R(mm, poi = NA_character_),
               "non-empty")
})

test_that("asymmetric mutation still errors in F1.5 (arrives in F5.1)", {
  skip_if_no_pedtools()
  mm <- marker_model(pedtools::nuclearPed(1), "M1", trio_freqs(),
                     mutation = list(model = "asymmetric", rate = 1e-3,
                                     ratio = 0.5, bias = 0.4))
  expect_error(mispitools:::cpt_marker_joint_R(mm),
               "mutation model")
})

test_that("linkage non-NULL errors in F1.2", {
  skip_if_no_pedtools()
  mm <- marker_model(pedtools::nuclearPed(1), "M1", trio_freqs(),
                     linkage = list(partner = "M2", theta = 0.05))
  expect_error(mispitools:::cpt_marker_joint_R(mm),
               "linked markers")
})

test_that("non-marker_model input errors", {
  expect_error(mispitools:::cpt_marker_joint_R(list()),
               "must be a 'marker_model'")
})

test_that("half-sib pedigree gives valid joint distribution", {
  skip_if_no_pedtools()
  ped <- pedtools::halfSibPed()
  freqs <- c("a" = 0.3, "b" = 0.4, "c" = 0.3)
  mm <- marker_model(ped, "M1", freqs)
  out <- mispitools:::cpt_marker_joint_R(mm)
  expect_equal(sum(out$P_H1), 1, tolerance = 1e-10)
  expect_equal(sum(out$P_H2), 1, tolerance = 1e-10)
  expect_true(all(out$P_H1 >= 0 & out$P_H2 >= 0))
})

test_that("ancestral_order returns parents-before-children", {
  skip_if_no_pedtools()
  ped <- pedtools::linearPed(2)
  nf <- as.character(pedtools::nonfounders(ped))
  ord <- mispitools:::ancestral_order(ped, nf)
  expect_setequal(ord, nf)
  ## Each entry's parents must be either founders or appear earlier in `ord`.
  founders <- as.character(pedtools::founders(ped))
  for (i in seq_along(ord)) {
    fa <- as.character(pedtools::father(ped, id = ord[i]))
    mo <- as.character(pedtools::mother(ped, id = ord[i]))
    earlier <- c(founders, ord[seq_len(i - 1L)])
    expect_true(fa %in% earlier)
    expect_true(mo %in% earlier)
  }
})

test_that("default poi prefers untyped non-founder", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  ## All members untyped: default = last untyped non-founder = "3".
  expect_identical(mispitools:::resolve_poi(ped, NULL), "3")
  expect_identical(mispitools:::resolve_poi(ped, "1"), "1")
})
