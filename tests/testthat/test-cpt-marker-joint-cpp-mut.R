## F2.4 — cpt_marker_joint_cpp() with mutation = "equal" / "stepwise".
##
## Cross-checks the C++ kernel against the R-reference engine (F1.4,
## F1.5) bit-for-bit (1e-12 tol) across the pedigrees used in F2.2 and
## a range of mutation rates / step ratios. Also runs a few sanity
## checks (sums-to-1, support lift under equal-rate, continuity as
## rate -> 0) to keep the C++ side honest independent of the R-ref.

skip_if_no_pedtools <- function() {
  testthat::skip_if_not_installed("pedtools")
}

compare_joint_mut <- function(model, tol = 1e-12, info = NULL) {
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
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Equal-rate cross-check vs R-ref (1e-12)
# ---------------------------------------------------------------------------

test_that("trio K=2 equal R=0.01 matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.6)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                     mutation = list(model = "equal", rate = 0.01))
  compare_joint_mut(mm, info = "trio K=2 equal R=0.01")
})

test_that("trio K=3 equal R=0.005 matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.3, "b" = 0.5, "c" = 0.2)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                     mutation = list(model = "equal", rate = 0.005))
  compare_joint_mut(mm, info = "trio K=3 equal R=0.005")
})

test_that("trio K=4 equal R=0.02 matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.1, "b" = 0.2, "c" = 0.3, "d" = 0.4)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                     mutation = list(model = "equal", rate = 0.02))
  compare_joint_mut(mm, info = "trio K=4 equal R=0.02")
})

test_that("halfSibPed K=3 equal R=0.01 matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.3, "c" = 0.3)
  mm <- marker_model(pedtools::halfSibPed(), "M1", freqs,
                     mutation = list(model = "equal", rate = 0.01))
  compare_joint_mut(mm, info = "halfSibPed K=3 equal R=0.01")
})

test_that("linearPed(2) K=2 equal R=0.005 matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.5, "b" = 0.5)
  mm <- marker_model(pedtools::linearPed(2), "M1", freqs,
                     mutation = list(model = "equal", rate = 0.005))
  compare_joint_mut(mm, info = "linearPed(2) K=2 equal R=0.005")
})

test_that("linearPed(2) K=3 equal R=0.02 matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.3, "b" = 0.4, "c" = 0.3)
  mm <- marker_model(pedtools::linearPed(2), "M1", freqs,
                     mutation = list(model = "equal", rate = 0.02))
  compare_joint_mut(mm, info = "linearPed(2) K=3 equal R=0.02")
})

test_that("nuclearPed(2) K=2 equal R=0.01 matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.6)
  mm <- marker_model(pedtools::nuclearPed(2), "M1", freqs,
                     mutation = list(model = "equal", rate = 0.01))
  compare_joint_mut(mm, info = "nuclearPed(2) K=2 equal R=0.01")
})

test_that("equal mutation P_H1 / P_H2 each sum to 1 (trio K=3, R=0.005)", {
  skip_if_no_pedtools()
  mm <- marker_model(pedtools::nuclearPed(1), "M1",
                     c("a" = 0.3, "b" = 0.5, "c" = 0.2),
                     mutation = list(model = "equal", rate = 0.005))
  out <- mispitools:::cpt_marker_joint_cpp_wrap(mm)
  expect_equal(sum(out$P_H1), 1, tolerance = 1e-12)
  expect_equal(sum(out$P_H2), 1, tolerance = 1e-12)
})

test_that("equal mutation lifts Mendelian zeros on AA x AB child (trio K=2)", {
  skip_if_no_pedtools()
  ## Under equal R>0, the (AA, AB) -> BB child is no longer impossible.
  freqs <- c("a" = 0.4, "b" = 0.6)
  R <- 0.05
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                     mutation = list(model = "equal", rate = R))
  out <- mispitools:::cpt_marker_joint_cpp_wrap(mm)
  rBB <- out[out$`1` == "a/a" & out$`2` == "a/b" & out$`3` == "b/b", ]
  expect_equal(nrow(rBB), 1L)
  expect_gt(rBB$P_H1, 0)
})

# ---------------------------------------------------------------------------
# Stepwise cross-check vs R-ref (1e-12)
# ---------------------------------------------------------------------------

test_that("trio K=2 stepwise (R=0.005, r=0.1) matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("12" = 0.4, "13" = 0.6)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                     mutation = list(model = "stepwise",
                                     rate = 0.005, ratio = 0.1))
  compare_joint_mut(mm, info = "trio K=2 stepwise")
})

test_that("trio K=3 stepwise unit-spaced (R=0.01, r=0.5) matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("10" = 0.3, "11" = 0.4, "12" = 0.3)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                     mutation = list(model = "stepwise",
                                     rate = 0.01, ratio = 0.5))
  compare_joint_mut(mm, info = "trio K=3 stepwise unit-spaced")
})

test_that("trio K=3 stepwise fractional labels (R=0.01, r=0.3)", {
  skip_if_no_pedtools()
  ## Microvariant-like fractional repeat counts; exercises the
  ## numeric_labels plumbing through the binding.
  freqs <- c("15" = 0.3, "15.2" = 0.3, "16" = 0.4)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                     mutation = list(model = "stepwise",
                                     rate = 0.01, ratio = 0.3))
  compare_joint_mut(mm, info = "trio K=3 fractional")
})

test_that("halfSibPed K=3 stepwise (R=0.005, r=0.1) matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("10" = 0.4, "11" = 0.3, "12" = 0.3)
  mm <- marker_model(pedtools::halfSibPed(), "M1", freqs,
                     mutation = list(model = "stepwise",
                                     rate = 0.005, ratio = 0.1))
  compare_joint_mut(mm, info = "halfSibPed K=3 stepwise")
})

test_that("linearPed(2) K=3 stepwise (R=0.02, r=0.2) matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("9" = 0.3, "10" = 0.4, "11" = 0.3)
  mm <- marker_model(pedtools::linearPed(2), "M1", freqs,
                     mutation = list(model = "stepwise",
                                     rate = 0.02, ratio = 0.2))
  compare_joint_mut(mm, info = "linearPed(2) K=3 stepwise")
})

test_that("nuclearPed(2) K=2 stepwise (R=0.005, r=0.1) matches R-ref to 1e-12", {
  skip_if_no_pedtools()
  freqs <- c("12" = 0.4, "13" = 0.6)
  mm <- marker_model(pedtools::nuclearPed(2), "M1", freqs,
                     mutation = list(model = "stepwise",
                                     rate = 0.005, ratio = 0.1))
  compare_joint_mut(mm, info = "nuclearPed(2) K=2 stepwise")
})

test_that("stepwise mutation P_H1 / P_H2 each sum to 1 (halfSib K=3)", {
  skip_if_no_pedtools()
  freqs <- c("10" = 0.4, "11" = 0.3, "12" = 0.3)
  mm <- marker_model(pedtools::halfSibPed(), "M1", freqs,
                     mutation = list(model = "stepwise",
                                     rate = 0.005, ratio = 0.1))
  out <- mispitools:::cpt_marker_joint_cpp_wrap(mm)
  expect_equal(sum(out$P_H1), 1, tolerance = 1e-12)
  expect_equal(sum(out$P_H2), 1, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Equal/stepwise vs none continuity and non-numeric label rejection
# ---------------------------------------------------------------------------

test_that("equal R -> 0 converges to mutation=none (trio K=3)", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.3, "b" = 0.5, "c" = 0.2)
  ped <- pedtools::nuclearPed(1)
  mm_none <- marker_model(ped, "M1", freqs)
  mm_eq   <- marker_model(ped, "M1", freqs,
                          mutation = list(model = "equal", rate = 1e-10))
  out_n <- mispitools:::cpt_marker_joint_cpp_wrap(mm_none)
  out_e <- mispitools:::cpt_marker_joint_cpp_wrap(mm_eq)

  ids <- setdiff(names(out_n), c("P_H1", "P_H2"))
  key_n <- do.call(paste, c(out_n[, ids, drop = FALSE], list(sep = "|")))
  key_e <- do.call(paste, c(out_e[, ids, drop = FALSE], list(sep = "|")))
  pn <- stats::setNames(out_n$P_H1, key_n)
  pe <- stats::setNames(out_e$P_H1, key_e)
  all_keys <- union(names(pn), names(pe))
  vn <- unname(pn[all_keys]); vn[is.na(vn)] <- 0
  ve <- unname(pe[all_keys]); ve[is.na(ve)] <- 0
  expect_lt(max(abs(vn - ve)), 1e-7)
})

test_that("stepwise wrap rejects non-numeric allele labels", {
  skip_if_no_pedtools()
  freqs <- c("A" = 0.4, "B" = 0.6)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                     mutation = list(model = "stepwise",
                                     rate = 1e-3, ratio = 0.1))
  expect_error(mispitools:::cpt_marker_joint_cpp_wrap(mm),
               "numeric allele labels")
})
