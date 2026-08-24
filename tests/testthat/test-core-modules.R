## Smoke tests for modules 1-4 (non-genetic LR, simulation, decision, combination).

# --- Module 1: Non-genetic LR functions ---

test_that("lr_sex returns a data.frame with expected columns", {
  res <- lr_sex(MPs = "M", eps = 0.05, numsims = 20, seed = 1)
  expect_s3_class(res, "data.frame")
  expect_true("Sexo" %in% names(res))
  expect_equal(nrow(res), 20)
})

test_that("lr_sex with LR=TRUE includes LRs column", {
  res <- lr_sex(MPs = "F", eps = 0.05, numsims = 20, LR = TRUE, seed = 1)
  expect_true("LRs" %in% names(res))
  expect_true(all(res$LRs > 0))
})

test_that("lr_age returns a data.frame with expected rows", {
  res <- lr_age(MPa = 30, MPr = 5, numsims = 20, seed = 1)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 20)
})

test_that("lr_hair_color returns a data.frame", {
  res <- lr_hair_color(MPc = 1, numsims = 20, seed = 1)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 20)
})

# --- Module 2: Simulation wrappers ---

test_that("sim_lr_prelim('sex') returns Related/Unrelated columns", {
  res <- sim_lr_prelim("sex", numsims = 20, seed = 1)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("Unrelated", "Related") %in% names(res)))
  expect_equal(nrow(res), 20)
})

test_that("sim_lr_prelim('age') returns correct shape", {
  res <- sim_lr_prelim("age", numsims = 20, seed = 1)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("Unrelated", "Related") %in% names(res)))
})

test_that("sim_lr_prelim rejects invalid vartype", {
  expect_error(sim_lr_prelim("invalid_type"), "vartype must be")
})

test_that("sim_lr_genetic returns Related and Unrelated lists", {
  skip_on_cran()
  skip_if_not_installed("forrel")
  suppressPackageStartupMessages({
    library(pedtools)
    library(forrel)
  })
  x <- linearPed(2)
  x <- setMarkers(x, locusAttributes = NorwegianFrequencies[1:3])
  x <- profileSim(x, N = 1, ids = 2, seed = 1)

  res <- sim_lr_genetic(x, missing = 5, numsims = 5, seed = 42)
  expect_type(res, "list")
  expect_true(all(c("Unrelated", "Related") %in% names(res)))
  expect_equal(length(res$Related), 5)
  expect_equal(length(res$Unrelated), 5)
})

test_that("lr_to_dataframe produces Related/Unrelated data.frame", {
  skip_on_cran()
  skip_if_not_installed("forrel")
  suppressPackageStartupMessages({
    library(pedtools)
    library(forrel)
  })
  x <- linearPed(2)
  x <- setMarkers(x, locusAttributes = NorwegianFrequencies[1:3])
  x <- profileSim(x, N = 1, ids = 2, seed = 1)

  res <- sim_lr_genetic(x, missing = 5, numsims = 5, seed = 42)
  df <- lr_to_dataframe(res)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("Unrelated", "Related") %in% names(df)))
  expect_equal(nrow(df), 5)
  expect_true(all(is.numeric(df$Related)))
})

# --- Module 3: Decision thresholds ---

test_that("decision_threshold returns a numeric scalar", {
  set.seed(1)
  sim_df <- data.frame(
    Unrelated = 10^runif(50, -3, 0),
    Related   = 10^runif(50,  2, 6)
  )
  thr <- decision_threshold(sim_df, weight = 10)
  expect_true(is.numeric(thr))
  expect_length(thr, 1)
})

test_that("threshold_rates returns a list", {
  set.seed(1)
  sim_df <- data.frame(
    Unrelated = 10^runif(50, -3, 0),
    Related   = 10^runif(50,  2, 6)
  )
  thr <- decision_threshold(sim_df, weight = 10)
  rates <- threshold_rates(sim_df, thr)
  expect_type(rates, "list")
})

# --- Module 4: Evidence combination ---

test_that("lr_combine multiplies paired draws and preserves nrow", {
  df1 <- data.frame(Unrelated = c(0.1, 0.2), Related = c(10, 20))
  df2 <- data.frame(Unrelated = c(2, 3),     Related = c(5, 7))
  res <- lr_combine(df1, df2)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_true(all(c("Unrelated", "Related") %in% names(res)))
  expect_equal(res$Unrelated[1], 0.1 * 2)
  expect_equal(res$Related[2], 20 * 7)
})

test_that("lr_combine rejects mismatched nrow", {
  df1 <- data.frame(Unrelated = c(0.1, 0.2), Related = c(10, 20))
  df2 <- data.frame(Unrelated = c(2),         Related = c(5))
  expect_error(lr_combine(df1, df2), "same number of rows")
})

test_that("lr_combine rejects non-data.frame", {
  expect_error(lr_combine("not a df", data.frame(Unrelated=1, Related=1)), "data.frame")
})
