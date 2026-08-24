## F1.1 — marker_model() S3 constructor: input validation + print method.
##
## The constructor itself does no joint computation; these tests pin the
## structural contract that the F1.2+ engines depend on.

skip_if_no_pedtools <- function() {
  testthat::skip_if_not_installed("pedtools")
}

make_ped <- function() {
  skip_if_no_pedtools()
  pedtools::nuclearPed(1)
}

valid_freqs <- function() {
  c("12" = 0.2, "13" = 0.3, "14" = 0.5)
}

test_that("constructor returns a marker_model object with stored inputs", {
  ped <- make_ped()
  freqs <- valid_freqs()
  mm <- marker_model(ped, "M1", freqs)

  expect_s3_class(mm, "marker_model")
  expect_identical(mm$marker_id, "M1")
  expect_identical(mm$freqs, freqs)
  expect_identical(mm$alleles, names(freqs))
  expect_identical(mm$mutation, list(model = "none", rate = 0))
  expect_null(mm$linkage)
  expect_true(pedtools::is.ped(mm$ped))
})

test_that("constructor stores supplied mutation and linkage lists", {
  ped <- make_ped()
  mm <- marker_model(
    ped, "M1", valid_freqs(),
    mutation = list(model = "stepwise", rate = 1e-3, ratio = 0.1),
    linkage = list(partner = "M2", theta = 0.05)
  )
  expect_identical(mm$mutation$model, "stepwise")
  expect_identical(mm$mutation$rate, 1e-3)
  expect_identical(mm$mutation$ratio, 0.1)
  expect_identical(mm$linkage$partner, "M2")
  expect_identical(mm$linkage$theta, 0.05)
})

test_that("invalid `ped` errors", {
  expect_error(
    marker_model("not a ped", "M1", valid_freqs()),
    "must be a 'ped' object"
  )
  expect_error(
    marker_model(NULL, "M1", valid_freqs()),
    "must be a 'ped' object"
  )
})

test_that("invalid `marker_id` errors", {
  ped <- make_ped()
  expect_error(marker_model(ped, NA_character_, valid_freqs()),
               "non-NA character")
  expect_error(marker_model(ped, c("a", "b"), valid_freqs()),
               "non-NA character")
  expect_error(marker_model(ped, 1, valid_freqs()), "non-NA character")
  expect_error(marker_model(ped, "", valid_freqs()), "non-empty")
})

test_that("invalid `freqs` errors", {
  ped <- make_ped()

  expect_error(marker_model(ped, "M1", c(0.5, 0.5)), "must be named")
  expect_error(marker_model(ped, "M1", c(a = 1.0)), "at least two")
  expect_error(
    marker_model(ped, "M1", c(a = 0.4, a = 0.6)),
    "duplicated allele labels"
  )
  expect_error(
    marker_model(ped, "M1", c(a = 0.4, b = NA_real_)),
    "must not contain NA"
  )
  expect_error(
    marker_model(ped, "M1", c(a = -0.1, b = 1.1)),
    "lie in \\[0, 1\\]"
  )
  expect_error(
    marker_model(ped, "M1", c(a = 0.3, b = 0.3)),
    "must sum to 1"
  )
  expect_error(
    marker_model(ped, "M1", "not numeric"),
    "must be a numeric vector"
  )
})

test_that("constructor accepts freqs that sum within tol", {
  ped <- make_ped()
  freqs <- c("a" = 0.5 - 5e-7, "b" = 0.5 + 5e-7)
  expect_no_error(marker_model(ped, "M1", freqs))
})

test_that("invalid `mutation` errors", {
  ped <- make_ped()
  expect_error(marker_model(ped, "M1", valid_freqs(), mutation = list()),
               "named list")
  expect_error(marker_model(ped, "M1", valid_freqs(),
                            mutation = list(rate = 0)),
               "must contain a `model` field")
  expect_error(marker_model(ped, "M1", valid_freqs(),
                            mutation = list(model = "weird")),
               "must be one of")
  expect_error(marker_model(ped, "M1", valid_freqs(),
                            mutation = list(model = "none", rate = 0.1)),
               "must be 0 when model = \"none\"")
  expect_error(marker_model(ped, "M1", valid_freqs(),
                            mutation = list(model = "equal")),
               "is required when model = \"equal\"")
  expect_error(marker_model(ped, "M1", valid_freqs(),
                            mutation = list(model = "equal", rate = 1)),
               "in \\[0, 1\\)")
  expect_error(marker_model(ped, "M1", valid_freqs(),
                            mutation = list(model = "equal", rate = -0.1)),
               "in \\[0, 1\\)")
  expect_error(
    marker_model(ped, "M1", valid_freqs(),
                 mutation = list(model = "stepwise", rate = 1e-3, ratio = 1.5)),
    "in \\(0, 1\\)"
  )
  expect_error(
    marker_model(ped, "M1", valid_freqs(),
                 mutation = list(model = "asymmetric", rate = 1e-3, bias = 2)),
    "in \\[0, 1\\]"
  )
})

test_that("invalid `linkage` errors", {
  ped <- make_ped()
  expect_error(
    marker_model(ped, "M1", valid_freqs(), linkage = list()),
    "named list"
  )
  expect_error(
    marker_model(ped, "M1", valid_freqs(),
                 linkage = list(theta = 0.05)),
    "linkage\\$partner"
  )
  expect_error(
    marker_model(ped, "M1", valid_freqs(),
                 linkage = list(partner = "M2", theta = -0.1)),
    "linkage\\$theta"
  )
  expect_error(
    marker_model(ped, "M1", valid_freqs(),
                 linkage = list(partner = "M2", theta = 0.6)),
    "linkage\\$theta"
  )
  expect_error(
    marker_model(ped, "M1", valid_freqs(),
                 linkage = list(partner = "", theta = 0.05)),
    "linkage\\$partner"
  )
})

test_that("print.marker_model emits a one-line-per-field summary", {
  ped <- make_ped()
  mm <- marker_model(ped, "M1", valid_freqs(),
                     mutation = list(model = "equal", rate = 1e-3))
  out <- capture.output(print(mm))
  expect_true(any(grepl("<marker_model>", out, fixed = TRUE)))
  expect_true(any(grepl("marker_id : M1", out, fixed = TRUE)))
  expect_true(any(grepl("alleles   : 3", out, fixed = TRUE)))
  expect_true(any(grepl("mutation  : equal", out)))
  expect_true(any(grepl("linkage   : none", out)))
})

test_that("print.marker_model returns its input invisibly", {
  ped <- make_ped()
  mm <- marker_model(ped, "M1", valid_freqs())
  expect_identical(withVisible(print(mm))$value, mm)
})
