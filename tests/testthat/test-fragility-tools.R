## Tests for calibrate_concentration_cutoff() and fragility_report().

test_that("fragility_report returns the expected fields for a balanced case", {
  lrs <- setNames(rep(2, 15), paste0("M", 1:15))
  fr <- fragility_report(lrs)

  expect_named(fr, c("total_log10_lr", "concentration", "top_marker",
                     "top_log10_lr", "residual_log10_lr", "flag",
                     "statement"))
  expect_equal(fr$total_log10_lr, 15 * log10(2))
  expect_equal(fr$concentration, 1 / 15)
  expect_true(grepl("W = ", fr$statement))
  expect_true(is.na(fr$flag))
})

test_that("fragility_report concentrates on the dominant marker", {
  lrs <- setNames(c(1e6, rep(1.05, 14)), c("D3", paste0("M", 1:14)))
  fr <- fragility_report(lrs)

  expect_equal(fr$top_marker, "D3")
  expect_gt(fr$concentration, 0.9)
  expect_lt(fr$residual_log10_lr, 1)
})

test_that("fragility_report flags cases above the cutoff", {
  lrs <- setNames(c(1e6, rep(1.05, 14)), c("D3", paste0("M", 1:14)))
  fr <- fragility_report(lrs, cutoff = 0.20, probs = 0.90)

  expect_true(fr$flag)
  expect_true(grepl("EXCEEDS", fr$statement))
  expect_true(grepl("90th-percentile", fr$statement))
})

test_that("fragility_report does not flag balanced cases below the cutoff", {
  lrs <- setNames(rep(2, 15), paste0("M", 1:15))
  fr <- fragility_report(lrs, cutoff = 0.20, probs = 0.90)

  expect_false(fr$flag)
  expect_true(grepl("below", fr$statement))
})

test_that("fragility_report validates inputs", {
  expect_error(fragility_report("not numeric"), "numeric vector")
  expect_error(fragility_report(c(1, 0, 2)), "must be > 0")
  expect_error(fragility_report(c(1, 2), cutoff = "bad"), "single numeric")
})

test_that("fragility_report assigns default marker names when missing", {
  lrs <- c(2, 3, 5)
  fr <- fragility_report(lrs)
  expect_equal(fr$top_marker, "M3")
})

test_that("fragility_report residual equals (1 - C_W+) * W", {
  lrs <- setNames(c(10, 2, 3, 1.5, 1.2), paste0("M", 1:5))
  fr <- fragility_report(lrs)
  log_lrs <- log10(lrs)
  expected_residual <- (1 - max(log_lrs) / sum(log_lrs[log_lrs > 0])) * sum(log_lrs)
  expect_equal(fr$residual_log10_lr, expected_residual)
})

test_that("calibrate_concentration_cutoff returns valid output on a tiny sim", {
  skip_on_cran()
  skip_if_not_installed("pedtools")
  skip_if_not_installed("forrel")
  suppressPackageStartupMessages({
    library(pedtools)
    library(forrel)
  })

  x <- linearPed(2)
  x <- setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
  x <- profileSim(x, N = 1, ids = 2, seed = 1)

  cal <- calibrate_concentration_cutoff(
    x, missing = 5, numsims = 20, probs = 0.90, seed = 42
  )

  expect_named(cal, c("cutoff", "probs", "distribution",
                      "total_log10_lr", "numsims"))
  expect_true(is.numeric(cal$cutoff))
  expect_true(cal$cutoff >= 0 && cal$cutoff <= 1)
  expect_equal(cal$probs, 0.90)
  expect_true(length(cal$distribution) <= 20)
  expect_true(cal$numsims > 0)
})

test_that("calibrate_concentration_cutoff validates probs argument", {
  skip_if_not_installed("pedtools")
  x <- pedtools::linearPed(1)
  expect_error(
    calibrate_concentration_cutoff(x, missing = 3, probs = 1.2),
    "probs must be"
  )
  expect_error(
    calibrate_concentration_cutoff(x, missing = 3, probs = -0.1),
    "probs must be"
  )
})
