## F6.3 — nongenetic_cpt_cpp() vs ng_cpt_R() cross-check.
##
## Validates that the C++ kernel reproduces the R-reference non-genetic
## CPT bit-for-bit (1e-12 tol) across the three feature classes
## (categorical / continuous / date) and their sub-cases. The legacy
## self-oracle (lr_sex / lr_age / ...) is the verifier's job in F6.6;
## here the oracle is the R reference engine ng_cpt_R().

compare_ng_cpt <- function(feature, tol = 1e-12, info = NULL) {
  ref  <- mispitools:::ng_cpt_R(feature)
  test <- mispitools:::ng_cpt_cpp_wrap(feature)

  expect_equal(nrow(test), nrow(ref), info = info)
  expect_identical(test$state, ref$state, info = info)
  expect_equal(test$p_h1, ref$p_h1, tolerance = tol, info = info)
  expect_equal(test$p_h2, ref$p_h2, tolerance = tol, info = info)
  expect_identical(attr(test, "feature_type"),
                   attr(ref, "feature_type"), info = info)
  expect_identical(attr(test, "observed"),
                   attr(ref, "observed"), info = info)
  invisible(NULL)
}

test_that("categorical scalar-eps, marginal reference, matches R-ref", {
  f <- nongenetic_feature(
    type = "sex", observed = "F",
    db_or_freqs = c(F = 0.55, M = 0.45), error = 0.05)
  compare_ng_cpt(f, info = "sex scalar marginal")
})

test_that("categorical scalar-eps, uniform reference, matches R-ref", {
  f <- nongenetic_feature(
    type = "region", observed = "NW",
    model = list(reference = "uniform"),
    db_or_freqs = c(NW = 0.2, NE = 0.3, S = 0.5), error = 0.1)
  compare_ng_cpt(f, info = "region scalar uniform")
})

test_that("categorical confusion matrix (positional) matches R-ref", {
  E <- matrix(c(0.80, 0.10, 0.10,
                0.05, 0.90, 0.05,
                0.15, 0.15, 0.70),
              nrow = 3, byrow = TRUE)
  f <- nongenetic_feature(
    type = "hair", observed = 2,
    db_or_freqs = c("1" = 0.3, "2" = 0.45, "3" = 0.25), error = E)
  compare_ng_cpt(f, info = "hair confusion matrix")
})

test_that("categorical asymmetric confusion matrix, observed first cat", {
  E <- matrix(c(0.7, 0.3,
                0.2, 0.8),
              nrow = 2, byrow = TRUE)
  f <- nongenetic_feature(
    type = "pigmentation", observed = "light",
    db_or_freqs = c(light = 0.4, dark = 0.6), error = E)
  compare_ng_cpt(f, info = "pigmentation 2x2 matrix")
})

test_that("continuous uniform reference matches R-ref", {
  f <- nongenetic_feature(
    type = "age", observed = 42,
    model = list(reference = "uniform", range = c(1, 80)),
    error = 0.05)
  compare_ng_cpt(f, info = "age uniform [1,80]")
})

test_that("continuous uniform with fractional range matches R-ref", {
  f <- nongenetic_feature(
    type = "age", observed = 10.4,
    model = list(reference = "uniform", range = c(2.3, 17.8)),
    error = 0.2)
  compare_ng_cpt(f, info = "age uniform fractional range")
})

test_that("continuous empirical reference matches R-ref", {
  set.seed(20260517)
  s <- c(5, 5, 12, 12, 12, 30, 30, 41, 41, 41, 41, 60, 7, 7, 7)
  f <- nongenetic_feature(
    type = "age", observed = 11,
    model = list(reference = "empirical"),
    db_or_freqs = s, error = 0.08)
  compare_ng_cpt(f, info = "age empirical")
})

test_that("continuous empirical with a tie in nearest cell matches R-ref", {
  s <- c(10, 20, 30, 40)
  ## observed exactly midway between 20 and 30 -> which.min takes the
  ## first (20); the C++ engine must use the same tie rule.
  f <- nongenetic_feature(
    type = "age", observed = 25,
    model = list(reference = "empirical"),
    db_or_freqs = s, error = 0.1)
  compare_ng_cpt(f, info = "age empirical tie")
})

test_that("date open search matches R-ref", {
  f <- nongenetic_feature(
    type = "birthdate", observed = 45,
    error = c(1, 4, 60, 11, 6, 4, 4))
  compare_ng_cpt(f, info = "birthdate open default")
})

test_that("date closed search with pre-binned freqs matches R-ref", {
  cuts <- c(-120, -30, 30, 120, 240, 360)
  f <- nongenetic_feature(
    type = "birthdate", observed = -10,
    model = list(search = "closed", cuts = cuts),
    db_or_freqs = c(3, 5, 40, 12, 8, 5, 7),
    error = c(2, 3, 50, 10, 7, 5, 6))
  compare_ng_cpt(f, info = "birthdate closed pre-binned")
})

test_that("custom feature routed through its declared class matches R-ref", {
  f <- nongenetic_feature(
    type = "custom", observed = "b",
    model = list(class = "categorical", reference = "marginal"),
    db_or_freqs = c(a = 0.2, b = 0.3, c = 0.5), error = 0.07)
  compare_ng_cpt(f, info = "custom categorical")
})

test_that("nongenetic_cpt_cpp errors on a bad feature class", {
  expect_error(
    mispitools:::nongenetic_cpt_cpp(
      feature_class = 9L, n_categories = 2L,
      error_is_matrix = FALSE, error_matrix = numeric(0),
      error_scalar = 0.05, observed_index = 0L,
      reference_uniform = FALSE, reference_freqs = c(0.5, 0.5),
      range_lo = 0, range_hi = 0, sample = numeric(0),
      observed_value = 0, n_bins = 0L, alpha = numeric(0),
      search_open = TRUE),
    "feature_class out of range")
})
