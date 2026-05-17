## F2.2 — cpt_marker_joint_cpp() vs cpt_marker_joint_R() cross-check.
##
## Validates that the C++ kernel reproduces the R-reference joint CPT
## bit-for-bit (1e-12 tol) across the pedigrees and markers used in F1,
## under mutation = "none". Equal/Stepwise are deferred to F2.4 tests.

skip_if_no_pedtools <- function() {
  testthat::skip_if_not_installed("pedtools")
}

compare_joint <- function(model, tol = 1e-12, info = NULL) {
  ref  <- mispitools:::cpt_marker_joint_R(model)
  test <- mispitools:::cpt_marker_joint_cpp_wrap(model)

  expect_equal(nrow(test), nrow(ref), info = info)
  expect_equal(names(test), names(ref), info = info)

  member_cols <- setdiff(names(ref), c("P_H1", "P_H2"))
  for (m in member_cols) {
    expect_identical(test[[m]], ref[[m]],
                     info = paste0(info, " | member=", m))
  }
  expect_equal(test$P_H1, ref$P_H1, tolerance = tol, info = info)
  expect_equal(test$P_H2, ref$P_H2, tolerance = tol, info = info)

  expect_identical(attr(test, "marker_id"), attr(ref, "marker_id"))
  expect_identical(attr(test, "alleles"), attr(ref, "alleles"))
  expect_identical(attr(test, "poi"), attr(ref, "poi"))
  invisible(NULL)
}

test_that("trio K=2 mut=none matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.6)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs)
  compare_joint(mm, info = "trio K=2")
})

test_that("trio K=3 mut=none matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.3, "b" = 0.5, "c" = 0.2)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs)
  compare_joint(mm, info = "trio K=3")
})

test_that("trio K=4 mut=none matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.1, "b" = 0.2, "c" = 0.3, "d" = 0.4)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs)
  compare_joint(mm, info = "trio K=4")
})

test_that("halfSibPed mut=none matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.3, "c" = 0.3)
  mm <- marker_model(pedtools::halfSibPed(), "M1", freqs)
  compare_joint(mm, info = "halfSibPed K=3")
})

test_that("linearPed(2) mut=none matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.5, "b" = 0.5)
  mm <- marker_model(pedtools::linearPed(2), "M1", freqs)
  compare_joint(mm, info = "linearPed(2) K=2")
})

test_that("linearPed(2) K=3 mut=none matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.3, "b" = 0.4, "c" = 0.3)
  mm <- marker_model(pedtools::linearPed(2), "M1", freqs)
  compare_joint(mm, info = "linearPed(2) K=3")
})

test_that("nuclearPed(2) (two children) mut=none matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.6)
  mm <- marker_model(pedtools::nuclearPed(2), "M1", freqs)
  compare_joint(mm, info = "nuclearPed(2) K=2")
})

test_that("P_H1 and P_H2 from the C++ engine each sum to 1 over a trio", {
  skip_if_no_pedtools()
  mm <- marker_model(pedtools::nuclearPed(1), "M1",
                     c("a" = 0.4, "b" = 0.6))
  out <- mispitools:::cpt_marker_joint_cpp_wrap(mm)
  expect_equal(sum(out$P_H1), 1, tolerance = 1e-12)
  expect_equal(sum(out$P_H2), 1, tolerance = 1e-12)
})

test_that("AAxAB Mendelian no-mutation: child is AA or AB with prob 0.5", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.6)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs)
  out <- mispitools:::cpt_marker_joint_cpp_wrap(mm)
  pAA <- freqs[["a"]]^2
  pAB <- 2 * freqs[["a"]] * freqs[["b"]]
  rAA <- out[out$`1` == "a/a" & out$`2` == "a/b" & out$`3` == "a/a", ]
  rAB <- out[out$`1` == "a/a" & out$`2` == "a/b" & out$`3` == "a/b", ]
  rBB <- out[out$`1` == "a/a" & out$`2` == "a/b" & out$`3` == "b/b", ]
  expect_equal(rAA$P_H1, pAA * pAB * 0.5, tolerance = 1e-12)
  expect_equal(rAB$P_H1, pAA * pAB * 0.5, tolerance = 1e-12)
  expect_equal(rBB$P_H1, 0, tolerance = 1e-12)
  expect_gt(rBB$P_H2, 0)
})

test_that("cpt_marker_joint_cpp_wrap() does not yet route asymmetric (F5.1)", {
  skip_if_no_pedtools()
  mm <- marker_model(pedtools::nuclearPed(1), "M1",
                     c("12" = 0.5, "13" = 0.5),
                     mutation = list(model = "asymmetric", rate = 0.005,
                                     ratio = 0.1, bias = 0.5))
  expect_error(mispitools:::cpt_marker_joint_cpp_wrap(mm),
               "asymmetric|not yet routed")
})

test_that("cpt_marker_joint_cpp() raw binding accepts Asymmetric/Dawid (F5.1)", {
  skip_if_no_pedtools()
  res <- cpt_marker_joint_cpp(
    father = c(-1L, -1L, 0L),
    mother = c(-1L, -1L, 1L),
    poi = 2L,
    freqs = c(0.4, 0.6),
    mutation_kind = 4L,
    mutation_rate = 0.005,
    mutation_range = 0.1
  )
  expect_true(is.list(res))
  expect_true(all(c("P_H1", "P_H2") %in% names(res)))
  expect_equal(sum(res$P_H1), 1, tolerance = 1e-12)
  expect_equal(sum(res$P_H2), 1, tolerance = 1e-12)
  expect_true(all(res$P_H1 >= 0) && all(res$P_H2 >= 0))
})

test_that("cpt_marker_joint_cpp() raw binding rejects Proportional kind", {
  skip_if_no_pedtools()
  expect_error(
    cpt_marker_joint_cpp(
      father = c(-1L, -1L, 0L),
      mother = c(-1L, -1L, 1L),
      poi = 2L,
      freqs = c(0.5, 0.5),
      mutation_kind = 3L,
      mutation_rate = 0.005
    ),
    "Proportional"
  )
})
