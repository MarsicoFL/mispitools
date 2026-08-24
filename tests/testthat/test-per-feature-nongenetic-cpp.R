## F6.4 — nongenetic per-feature KL / LR-distribution C++ kernels vs the
## R reference engine (r_ref_nongenetic.R). The C++ path adapts the
## one-dimensional non-genetic CPT into a single-member JointTable and
## reuses the genetic per_marker_kl / per_marker_lr_dist kernels, so it
## must reproduce per_feature_kl_R / per_feature_lr_dist_R bit-for-bit.
## Tolerance 1e-12 (the F6.x cross-check spec). The legacy self-oracle
## (lr_sex / lr_age / ...) is the verifier's job in F6.6.

compare_pf_kl <- function(feature, tol = 1e-12, info = NULL) {
  ref  <- mispitools:::per_feature_kl_R(feature)
  test <- mispitools:::per_feature_kl_cpp_wrap(feature)

  expect_identical(test$feature, ref$feature, info = info)
  expect_equal(test$e_log10_lr_h1, ref$e_log10_lr_h1,
               tolerance = tol, info = info)
  expect_equal(test$e_log10_lr_h2, ref$e_log10_lr_h2,
               tolerance = tol, info = info)
  expect_equal(test$kl_h1_to_h2, ref$kl_h1_to_h2,
               tolerance = tol, info = info)
  expect_equal(test$kl_h2_to_h1, ref$kl_h2_to_h1,
               tolerance = tol, info = info)
  ## Self-consistency: kl_h1_to_h2 == e_log10_lr_h1 * ln(10).
  if (is.finite(test$e_log10_lr_h1)) {
    expect_equal(test$kl_h1_to_h2, test$e_log10_lr_h1 * log(10),
                 tolerance = tol, info = info)
  }
  invisible(NULL)
}

compare_pf_lr_dist <- function(feature, aggregate = TRUE,
                               tol = 1e-12, info = NULL) {
  ref  <- mispitools:::per_feature_lr_dist_R(feature, aggregate = aggregate)
  test <- mispitools:::per_feature_lr_dist_cpp_wrap(feature,
                                                    aggregate = aggregate)

  expect_equal(nrow(test), nrow(ref), info = info)
  expect_equal(test$log10_lr, ref$log10_lr, tolerance = tol, info = info)
  expect_equal(test$p_h1, ref$p_h1, tolerance = tol, info = info)
  expect_equal(test$p_h2, ref$p_h2, tolerance = tol, info = info)
  expect_identical(attr(test, "feature_type"),
                   attr(ref, "feature_type"), info = info)
  expect_identical(attr(test, "observed"),
                   attr(ref, "observed"), info = info)
  invisible(NULL)
}

## --- categorical -----------------------------------------------------

test_that("categorical scalar-eps marginal (sex): KL + LR dist", {
  f <- nongenetic_feature(
    type = "sex", observed = "F",
    db_or_freqs = c(F = 0.55, M = 0.45), error = 0.05)
  compare_pf_kl(f, info = "sex scalar marginal KL")
  compare_pf_lr_dist(f, info = "sex scalar marginal LR dist")
  compare_pf_lr_dist(f, aggregate = FALSE,
                     info = "sex scalar marginal LR dist raw")
})

test_that("categorical scalar-eps uniform reference (region)", {
  f <- nongenetic_feature(
    type = "region", observed = "NW",
    model = list(reference = "uniform"),
    db_or_freqs = c(NW = 0.2, NE = 0.3, S = 0.5), error = 0.1)
  compare_pf_kl(f, info = "region scalar uniform KL")
  compare_pf_lr_dist(f, info = "region scalar uniform LR dist")
})

test_that("categorical confusion matrix (hair)", {
  E <- matrix(c(0.80, 0.10, 0.10,
                0.05, 0.90, 0.05,
                0.15, 0.15, 0.70),
              nrow = 3, byrow = TRUE)
  f <- nongenetic_feature(
    type = "hair", observed = 2,
    db_or_freqs = c("1" = 0.3, "2" = 0.45, "3" = 0.25), error = E)
  compare_pf_kl(f, info = "hair confusion matrix KL")
  compare_pf_lr_dist(f, info = "hair confusion matrix LR dist")
  compare_pf_lr_dist(f, aggregate = FALSE,
                     info = "hair confusion matrix LR dist raw")
})

test_that("categorical asymmetric 2x2 matrix (pigmentation)", {
  E <- matrix(c(0.7, 0.3,
                0.2, 0.8),
              nrow = 2, byrow = TRUE)
  f <- nongenetic_feature(
    type = "pigmentation", observed = "light",
    db_or_freqs = c(light = 0.4, dark = 0.6), error = E)
  compare_pf_kl(f, info = "pigmentation 2x2 KL")
  compare_pf_lr_dist(f, info = "pigmentation 2x2 LR dist")
})

test_that("categorical with a zero-mass H1 cell propagates -Inf", {
  ## A deterministic confusion matrix (no error off the diagonal) puts
  ## zero H1 mass on the non-observed categories; with a uniform H2 each
  ## is a (p_h1 = 0, p_h2 > 0) atom -> log10 LR = -Inf, matching R-ref.
  E <- matrix(c(1, 0, 0,
                0, 1, 0,
                0, 0, 1),
              nrow = 3, byrow = TRUE)
  f <- nongenetic_feature(
    type = "hair", observed = 1,
    model = list(reference = "uniform"),
    db_or_freqs = c("1" = 0.5, "2" = 0.3, "3" = 0.2), error = E)
  ref  <- mispitools:::per_feature_kl_R(f)
  test <- mispitools:::per_feature_kl_cpp_wrap(f)
  expect_equal(test$e_log10_lr_h1, ref$e_log10_lr_h1)
  expect_true(is.infinite(test$e_log10_lr_h2) && test$e_log10_lr_h2 < 0)
  expect_identical(test$e_log10_lr_h2, ref$e_log10_lr_h2)
  expect_identical(test$kl_h2_to_h1, ref$kl_h2_to_h1)  # +Inf both

  d <- mispitools:::per_feature_lr_dist_cpp_wrap(f)
  dr <- mispitools:::per_feature_lr_dist_R(f)
  expect_equal(nrow(d), nrow(dr))
  expect_equal(d$log10_lr, dr$log10_lr)
  expect_equal(d$p_h1, dr$p_h1)
  expect_equal(d$p_h2, dr$p_h2)
})

## --- continuous ------------------------------------------------------

test_that("continuous uniform reference (age [1,80])", {
  f <- nongenetic_feature(
    type = "age", observed = 42,
    model = list(reference = "uniform", range = c(1, 80)),
    error = 0.05)
  compare_pf_kl(f, info = "age uniform KL")
  compare_pf_lr_dist(f, info = "age uniform LR dist")
  compare_pf_lr_dist(f, aggregate = FALSE,
                     info = "age uniform LR dist raw")
})

test_that("continuous empirical reference (age)", {
  s <- c(5, 5, 12, 12, 12, 30, 30, 41, 41, 41, 41, 60, 7, 7, 7)
  f <- nongenetic_feature(
    type = "age", observed = 11,
    model = list(reference = "empirical"),
    db_or_freqs = s, error = 0.08)
  compare_pf_kl(f, info = "age empirical KL")
  compare_pf_lr_dist(f, info = "age empirical LR dist")
})

## --- date ------------------------------------------------------------

test_that("date open search (birthdate)", {
  f <- nongenetic_feature(
    type = "birthdate", observed = 45,
    error = c(1, 4, 60, 11, 6, 4, 4))
  compare_pf_kl(f, info = "birthdate open KL")
  compare_pf_lr_dist(f, info = "birthdate open LR dist")
})

test_that("date closed search with pre-binned freqs (birthdate)", {
  cuts <- c(-120, -30, 30, 120, 240, 360)
  f <- nongenetic_feature(
    type = "birthdate", observed = -10,
    model = list(search = "closed", cuts = cuts),
    db_or_freqs = c(3, 5, 40, 12, 8, 5, 7),
    error = c(2, 3, 50, 10, 7, 5, 6))
  compare_pf_kl(f, info = "birthdate closed KL")
  compare_pf_lr_dist(f, info = "birthdate closed LR dist")
})

test_that("custom feature routed through its declared class", {
  f <- nongenetic_feature(
    type = "custom", observed = "b",
    model = list(class = "categorical", reference = "marginal"),
    db_or_freqs = c(a = 0.2, b = 0.3, c = 0.5), error = 0.07)
  compare_pf_kl(f, info = "custom categorical KL")
  compare_pf_lr_dist(f, info = "custom categorical LR dist")
})

## --- error path ------------------------------------------------------

test_that("non-genetic per-feature kernels reject a bad feature class", {
  expect_error(
    mispitools:::nongenetic_per_feature_kl_cpp(
      feature_class = 9L, n_categories = 2L,
      error_is_matrix = FALSE, error_matrix = numeric(0),
      error_scalar = 0.05, observed_index = 0L,
      reference_uniform = FALSE, reference_freqs = c(0.5, 0.5),
      range_lo = 0, range_hi = 0, sample = numeric(0),
      observed_value = 0, n_bins = 0L, alpha = numeric(0),
      search_open = TRUE),
    "feature_class out of range")
  expect_error(
    mispitools:::nongenetic_per_feature_lr_dist_cpp(
      feature_class = 9L, n_categories = 2L,
      error_is_matrix = FALSE, error_matrix = numeric(0),
      error_scalar = 0.05, observed_index = 0L,
      reference_uniform = FALSE, reference_freqs = c(0.5, 0.5),
      range_lo = 0, range_hi = 0, sample = numeric(0),
      observed_value = 0, n_bins = 0L, alpha = numeric(0),
      search_open = TRUE, aggregate = TRUE),
    "feature_class out of range")
  ## A CPT construction error (too few categories) propagates from the
  ## kernel as a stop(), not a crash across the boundary.
  expect_error(
    mispitools:::nongenetic_per_feature_kl_cpp(
      feature_class = 0L, n_categories = 1L,
      error_is_matrix = FALSE, error_matrix = numeric(0),
      error_scalar = 0.05, observed_index = 0L,
      reference_uniform = TRUE, reference_freqs = numeric(0),
      range_lo = 0, range_hi = 0, sample = numeric(0),
      observed_value = 0, n_bins = 0L, alpha = numeric(0),
      search_open = TRUE),
    "at least two categories")
})
