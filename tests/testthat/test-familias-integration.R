## Tests for Familias integration wrapper

make_fake_familias <- function(lrs, ref_name = "Unrelated",
                                test_name = "RelatedPedigree",
                                marker_names = NULL) {
  if (is.null(marker_names)) marker_names <- paste0("M", seq_along(lrs))
  mat <- cbind(rep(1, length(lrs)), unname(lrs))
  dimnames(mat) <- list(marker_names, c(ref_name, test_name))
  list(LRperMarker = mat)
}

# ---- Structural checks ----

test_that("familias_trajectory returns a list with expected components", {
  fam <- make_fake_familias(c(5, 2, 8, 3, 1.5, 0.8))
  out <- familias_trajectory(fam, test_pedigree = 2)
  expect_type(out, "list")
  expect_named(out, c("lrs", "trajectory", "metrics",
                       "concentration_positive", "leave_one_out"))
  expect_s3_class(out$trajectory, "data.frame")
  expect_s3_class(out$leave_one_out, "data.frame")
  expect_type(out$metrics, "list")
})

test_that("familias_trajectory preserves marker names", {
  lrs <- c(D3S1358 = 5, TH01 = 2, D21S11 = 8)
  fam <- make_fake_familias(lrs, marker_names = names(lrs))
  out <- familias_trajectory(fam, test_pedigree = 2)
  expect_equal(names(out$lrs), names(lrs))
  expect_equal(unname(out$lrs), unname(lrs))
})

test_that("familias_trajectory accepts test_pedigree by name", {
  fam <- make_fake_familias(c(5, 2, 8), test_name = "UncleNephew")
  out <- familias_trajectory(fam, test_pedigree = "UncleNephew")
  expect_equal(length(out$lrs), 3)
})

# ---- Consistency with underlying primitives ----

test_that("concentration_positive equals concentration_index_positive on log10(LR)", {
  lrs <- c(5, 2, 80, 3, 1.5)
  fam <- make_fake_familias(lrs)
  out <- familias_trajectory(fam, test_pedigree = 2)
  expect_equal(out$concentration_positive,
               concentration_index_positive(log10(lrs)),
               tolerance = 1e-12)
})

test_that("trajectory posterior at final step matches cumulative LR / (1 + cumulative LR)", {
  lrs <- c(5, 2, 8)
  fam <- make_fake_familias(lrs)
  out <- familias_trajectory(fam, test_pedigree = 2)
  cum_lr <- prod(lrs)
  expected_posterior_h1 <- cum_lr / (1 + cum_lr)
  expect_equal(tail(out$trajectory$posterior_h1, 1),
               expected_posterior_h1,
               tolerance = 1e-12)
})

# ---- Error handling ----

test_that("familias_trajectory errors on invalid input", {
  expect_error(familias_trajectory("not a list"), "must be a list")
  expect_error(familias_trajectory(list()), "must contain a matrix named 'LRperMarker'")
  expect_error(familias_trajectory(list(LRperMarker = 1:5)), "numeric matrix")
  expect_error(
    familias_trajectory(list(LRperMarker = matrix(1, nrow = 3, ncol = 1))),
    "at least 2 columns"
  )
})

test_that("familias_trajectory errors on non-positive LRs", {
  fam <- make_fake_familias(c(5, -1, 3))
  expect_error(familias_trajectory(fam, test_pedigree = 2), "strictly positive")
})

test_that("familias_trajectory errors when test_pedigree equals ref_pedigree", {
  fam <- make_fake_familias(c(5, 2, 8))
  expect_error(
    familias_trajectory(fam, test_pedigree = 1, ref_pedigree = 1),
    "must be different"
  )
})

test_that("familias_trajectory errors on unknown pedigree name", {
  fam <- make_fake_familias(c(5, 2, 8), test_name = "Related")
  expect_error(
    familias_trajectory(fam, test_pedigree = "Nonexistent"),
    "not found"
  )
})
