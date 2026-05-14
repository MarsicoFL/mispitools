## F3.2 — public R API: per_marker_kl() / per_marker_kl_profile().
##
## Thin R layer over cpt_marker_joint_cpp_wrap() + cpp_per_marker_kl().
## Tests pin the schema, cross-check against the F1.6 R reference engine
## bit-for-bit on finite cases, and exercise the profile loop / POI
## resolution / list-name overrides.

skip_if_no_pedtools <- function() {
  testthat::skip_if_not_installed("pedtools")
}

trio_model <- function(freqs = c("a" = 0.4, "b" = 0.6),
                       mutation = list(model = "none", rate = 0),
                       id = "M1") {
  skip_if_no_pedtools()
  marker_model(pedtools::nuclearPed(1), id, freqs, mutation = mutation)
}

expected_cols <- c("marker", "e_log10_lr_h1", "e_log10_lr_h2",
                   "kl_h1_to_h2", "kl_h2_to_h1",
                   "abs_cont_violations_h1", "abs_cont_violations_h2",
                   "mass_violations_h1", "mass_violations_h2")

# ---------------------------------------------------------------------------
# Structure
# ---------------------------------------------------------------------------

test_that("per_marker_kl returns the documented one-row schema", {
  mm <- trio_model(mutation = list(model = "equal", rate = 0.005))
  out <- per_marker_kl(mm)

  expect_s3_class(out, "data.frame")
  expect_named(out, expected_cols)
  expect_equal(nrow(out), 1L)
  expect_identical(out$marker, "M1")
  expect_true(is.numeric(out$e_log10_lr_h1))
  expect_true(is.numeric(out$kl_h1_to_h2))
  expect_true(is.integer(out$abs_cont_violations_h1))
  expect_identical(attr(out, "poi"), "3")
})

test_that("per_marker_kl rejects non-marker_model input", {
  expect_error(per_marker_kl(list(a = 1)), "must be a 'marker_model'")
  expect_error(per_marker_kl(NULL), "must be a 'marker_model'")
})

# ---------------------------------------------------------------------------
# Cross-check vs R-reference (finite cases bit-for-bit, ±Inf via identical)
# ---------------------------------------------------------------------------

expect_matches_ref <- function(out, ref, tol = 1e-12, info = NULL) {
  expect_identical(out$marker, ref$marker, info = info)
  expect_equal(out$e_log10_lr_h1, ref$e_log10_lr_h1, tolerance = tol, info = info)
  expect_equal(out$e_log10_lr_h2, ref$e_log10_lr_h2, tolerance = tol, info = info)
  expect_equal(out$kl_h1_to_h2,   ref$kl_h1_to_h2,   tolerance = tol, info = info)
  expect_equal(out$kl_h2_to_h1,   ref$kl_h2_to_h1,   tolerance = tol, info = info)
}

test_that("per_marker_kl matches R-ref on trio K=2 mut=equal", {
  skip_if_no_pedtools()
  mm <- trio_model(mutation = list(model = "equal", rate = 0.005))
  ref <- mispitools:::per_marker_kl_R(mm)
  out <- per_marker_kl(mm)
  expect_matches_ref(out, ref, info = "trio K=2 mut=equal")
})

test_that("per_marker_kl matches R-ref on trio K=3 mut=stepwise", {
  skip_if_no_pedtools()
  freqs <- c("12" = 0.3, "13" = 0.4, "14" = 0.3)
  mm <- marker_model(pedtools::nuclearPed(1), "M1", freqs,
                     mutation = list(model = "stepwise", rate = 0.005,
                                     ratio = 0.1))
  ref <- mispitools:::per_marker_kl_R(mm)
  out <- per_marker_kl(mm)
  expect_matches_ref(out, ref, info = "trio K=3 mut=stepwise")
  expect_equal(out$abs_cont_violations_h1, 0L)
  expect_equal(out$abs_cont_violations_h2, 0L)
})

test_that("per_marker_kl matches R-ref on halfSibPed mut=equal", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.3, "c" = 0.3)
  mm <- marker_model(pedtools::halfSibPed(), "MhalfSib", freqs,
                     mutation = list(model = "equal", rate = 0.01))
  ref <- mispitools:::per_marker_kl_R(mm)
  out <- per_marker_kl(mm)
  expect_matches_ref(out, ref, info = "halfSibPed K=3 mut=equal")
})

test_that("per_marker_kl matches R-ref on nuclearPed(2) mut=none (Inf branch)", {
  skip_if_no_pedtools()
  freqs <- c("a" = 0.4, "b" = 0.6)
  mm <- marker_model(pedtools::nuclearPed(2), "M1", freqs)
  ref <- mispitools:::per_marker_kl_R(mm)
  out <- per_marker_kl(mm)
  expect_equal(out$kl_h1_to_h2, ref$kl_h1_to_h2, tolerance = 1e-12)
  expect_identical(out$kl_h2_to_h1, ref$kl_h2_to_h1)
  expect_equal(out$e_log10_lr_h1, ref$e_log10_lr_h1, tolerance = 1e-12)
  expect_identical(out$e_log10_lr_h2, ref$e_log10_lr_h2)
  expect_gte(out$abs_cont_violations_h1, 1L)
  expect_equal(out$abs_cont_violations_h2, 0L)
  expect_gt(out$mass_violations_h1, 0)
})

test_that("per_marker_kl: trio HWE K=2 p=0.5 mut=none analytic", {
  skip_if_no_pedtools()
  mm <- trio_model(freqs = c("a" = 0.5, "b" = 0.5))
  out <- per_marker_kl(mm)
  expect_equal(out$e_log10_lr_h1, 0.625 * log10(2), tolerance = 1e-12)
  expect_equal(out$kl_h1_to_h2, 0.625 * log10(2) * log(10), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Self-consistency (KL = E[log_e LR] = E[log10 LR] * ln(10))
# ---------------------------------------------------------------------------

test_that("kl identities hold under mutation = equal", {
  mm <- trio_model(freqs = c("a" = 0.3, "b" = 0.7),
                   mutation = list(model = "equal", rate = 0.005))
  out <- per_marker_kl(mm)
  expect_equal(out$kl_h1_to_h2, out$e_log10_lr_h1 * log(10), tolerance = 1e-12)
  expect_equal(out$kl_h2_to_h1, -out$e_log10_lr_h2 * log(10), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# POI passthrough
# ---------------------------------------------------------------------------

test_that("per_marker_kl: explicit poi propagates and validates", {
  skip_if_no_pedtools()
  mm <- trio_model(mutation = list(model = "equal", rate = 0.005))
  out <- per_marker_kl(mm, poi = "3")
  expect_identical(attr(out, "poi"), "3")
  expect_error(per_marker_kl(mm, poi = "nonexistent_member"),
               "is not a member of the pedigree")
})

# ---------------------------------------------------------------------------
# per_marker_kl_profile
# ---------------------------------------------------------------------------

test_that("per_marker_kl_profile returns one row per model in order", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  models <- list(
    marker_model(ped, "M1", c("a" = 0.4, "b" = 0.6),
                 mutation = list(model = "equal", rate = 0.005)),
    marker_model(ped, "M2", c("a" = 0.2, "b" = 0.3, "c" = 0.5),
                 mutation = list(model = "equal", rate = 0.005))
  )
  out <- per_marker_kl_profile(models)
  expect_s3_class(out, "data.frame")
  expect_named(out, expected_cols)
  expect_equal(nrow(out), 2L)
  expect_identical(out$marker, c("M1", "M2"))
  expect_identical(rownames(out), c("1", "2"))
  expect_identical(attr(out, "poi"), c("3", "3"))
})

test_that("per_marker_kl_profile: rows agree with per_marker_kl(model)", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  models <- list(
    marker_model(ped, "M1", c("a" = 0.4, "b" = 0.6),
                 mutation = list(model = "equal", rate = 0.005)),
    marker_model(ped, "M2", c("12" = 0.2, "13" = 0.3, "14" = 0.5),
                 mutation = list(model = "stepwise", rate = 0.005,
                                 ratio = 0.1))
  )
  out <- per_marker_kl_profile(models)
  for (i in seq_along(models)) {
    s <- per_marker_kl(models[[i]])
    for (col in setdiff(expected_cols, "marker")) {
      expect_equal(out[[col]][i], s[[col]],
                   tolerance = 1e-15, info = paste(col, i))
    }
  }
})

test_that("per_marker_kl_profile: list names override marker_id only when non-empty", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  m1 <- marker_model(ped, "internal_A", c("a" = 0.4, "b" = 0.6),
                     mutation = list(model = "equal", rate = 0.005))
  m2 <- marker_model(ped, "internal_B", c("a" = 0.4, "b" = 0.6),
                     mutation = list(model = "equal", rate = 0.005))
  models <- list(m1, m2)
  names(models) <- c("D3S1358", "")
  out <- per_marker_kl_profile(models)
  expect_identical(out$marker, c("D3S1358", "internal_B"))
})

test_that("per_marker_kl_profile: scalar poi applied to all models", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  models <- list(
    marker_model(ped, "M1", c("a" = 0.4, "b" = 0.6),
                 mutation = list(model = "equal", rate = 0.005)),
    marker_model(ped, "M2", c("a" = 0.4, "b" = 0.6),
                 mutation = list(model = "equal", rate = 0.005))
  )
  out <- per_marker_kl_profile(models, poi = "3")
  expect_identical(attr(out, "poi"), "3")
  expect_equal(nrow(out), 2L)
})

test_that("per_marker_kl_profile: validation errors", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  mm <- marker_model(ped, "M1", c("a" = 0.4, "b" = 0.6),
                     mutation = list(model = "equal", rate = 0.005))
  expect_error(per_marker_kl_profile("not_a_list"),
               "must be a list of 'marker_model'")
  expect_error(per_marker_kl_profile(list()),
               "must contain at least one")
  expect_error(per_marker_kl_profile(list(mm, list(a = 1))),
               "non-conforming entries at positions: 2")
  expect_error(per_marker_kl_profile(list(mm), poi = c("3", "4")),
               "single non-empty character string")
  expect_error(per_marker_kl_profile(list(mm), poi = ""),
               "single non-empty character string")
})
