## F6.1 -- nongenetic_feature() S3 constructor: input validation + print.
##
## Design milestone: the constructor does no CPT computation. These tests
## pin the structural contract that the F6.2+ migration and the F6.3 C++
## engine depend on (symmetry with marker_model / F1.1).

test_that("categorical sex feature stores validated inputs", {
  f <- nongenetic_feature(
    type = "sex",
    observed = "F",
    db_or_freqs = c(F = 0.5, M = 0.5),
    error = 0.05
  )
  expect_s3_class(f, "nongenetic_feature")
  expect_identical(f$type, "sex")
  expect_identical(f$feature_class, "categorical")
  expect_identical(f$observed, "F")
  expect_identical(f$categories, c("F", "M"))
  expect_identical(f$model, list(reference = "marginal"))
  expect_identical(f$error, 0.05)
})

test_that("numeric categorical observed is matched as character", {
  f <- nongenetic_feature(
    type = "hair",
    observed = 1,
    db_or_freqs = c("1" = .3, "2" = .2, "3" = .25, "4" = .15, "5" = .1),
    error = error_matrix_hair()
  )
  expect_identical(f$observed, "1")
  expect_true(is.matrix(f$error))
  expect_equal(dim(f$error), c(5L, 5L))
})

test_that("confusion matrix is positional; dimnames are informational", {
  # error_matrix_hair() is label-keyed (Black/Brown/...) but categories
  # here are "1".."5"; the matrix must be accepted positionally.
  f <- nongenetic_feature(
    type = "hair",
    observed = 1,
    db_or_freqs = c("1" = .3, "2" = .2, "3" = .25, "4" = .15, "5" = .1),
    error = error_matrix_hair()
  )
  expect_true(is.matrix(f$error))
  expect_equal(unname(dim(f$error)), c(5L, 5L))

  # Non-square / non-row-stochastic matrices are still rejected.
  expect_error(
    nongenetic_feature(
      type = "hair", observed = 1,
      db_or_freqs = c("1" = .3, "2" = .2, "3" = .25, "4" = .15, "5" = .1),
      error = matrix(0.5, 3, 3)
    ),
    "must be 5 x 5"
  )
  expect_error(
    nongenetic_feature(
      type = "sex", observed = "F",
      db_or_freqs = c(F = .5, M = .5),
      error = matrix(c(0.9, 0.2, 0.1, 0.7), 2, 2)
    ),
    "rows must each sum to 1"
  )
})

test_that("continuous age feature: uniform reference defaults + storage", {
  f <- nongenetic_feature(
    type = "age",
    observed = 42,
    error = 0.05
  )
  expect_identical(f$feature_class, "continuous")
  expect_identical(f$model$reference, "uniform")
  expect_identical(f$model$range, c(1, 80))
  expect_null(f$db_or_freqs)
  expect_identical(f$observed, 42)
})

test_that("continuous empirical reference requires a numeric sample", {
  f <- nongenetic_feature(
    type = "age",
    observed = 30,
    model = list(reference = "empirical"),
    db_or_freqs = c(20, 25, 33, 41, 58),
    error = 0.05
  )
  expect_identical(f$model$reference, "empirical")
  expect_identical(f$db_or_freqs, c(20, 25, 33, 41, 58))

  expect_error(
    nongenetic_feature(
      type = "age", observed = 30,
      model = list(reference = "empirical"),
      db_or_freqs = NULL, error = 0.05
    ),
    "finite numeric sample"
  )
})

test_that("uniform continuous reference rejects a supplied db", {
  expect_error(
    nongenetic_feature(
      type = "age", observed = 40,
      model = list(reference = "uniform", range = c(1, 80)),
      db_or_freqs = c(1, 2, 3), error = 0.05
    ),
    "must be NULL when the continuous reference is"
  )
})

test_that("date birthdate feature: open search defaults", {
  f <- nongenetic_feature(
    type = "birthdate",
    observed = 45,
    error = c(1, 4, 60, 11, 6, 4, 4)
  )
  expect_identical(f$feature_class, "date")
  expect_identical(f$model$search, "open")
  expect_identical(f$model$cuts, c(-120, -30, 30, 120, 240, 360))
  expect_null(f$db_or_freqs)
  expect_identical(f$observed, 45)
})

test_that("date observed accepts Date and YYYY-MM-DD string", {
  f1 <- nongenetic_feature("birthdate", observed = "1976-07-15",
                           error = c(1, 4, 60, 11, 6, 4, 4))
  expect_s3_class(f1$observed, "Date")
  f2 <- nongenetic_feature("birthdate", observed = as.Date("1976-07-15"),
                           error = c(1, 4, 60, 11, 6, 4, 4))
  expect_s3_class(f2$observed, "Date")
  expect_error(
    nongenetic_feature("birthdate", observed = "not-a-date",
                       error = c(1, 4, 60, 11, 6, 4, 4)),
    "must be a Date"
  )
})

test_that("closed date search needs bin freqs of length(cuts)+1", {
  f <- nongenetic_feature(
    type = "birthdate", observed = 10,
    model = list(search = "closed", cuts = c(-30, 30)),
    db_or_freqs = c(2, 10, 3),
    error = c(1, 1, 1)
  )
  expect_identical(f$db_or_freqs, c(2, 10, 3))
  expect_error(
    nongenetic_feature(
      type = "birthdate", observed = 10,
      model = list(search = "closed", cuts = c(-30, 30)),
      db_or_freqs = c(2, 10),
      error = c(1, 1, 1)
    ),
    "closed date search must be"
  )
})

test_that("Dirichlet alpha length must match the bin count", {
  expect_error(
    nongenetic_feature(
      type = "birthdate", observed = 10,
      model = list(search = "open", cuts = c(-30, 30)),
      error = c(1, 1)
    ),
    "Dirichlet alpha vector of length 3"
  )
})

test_that("custom feature requires model$class and follows that class", {
  f <- nongenetic_feature(
    type = "custom",
    observed = "x",
    model = list(class = "categorical"),
    db_or_freqs = c(x = 0.4, y = 0.6),
    error = 0.1
  )
  expect_identical(f$feature_class, "categorical")
  expect_identical(f$categories, c("x", "y"))

  expect_error(
    nongenetic_feature(type = "custom", observed = "x",
                       db_or_freqs = c(x = 1), error = 0.1),
    "must be supplied"
  )
  expect_error(
    nongenetic_feature(type = "custom", observed = "x",
                       model = list(class = "bogus"),
                       db_or_freqs = c(x = .5, y = .5), error = 0.1),
    "model\\$class` must be one of"
  )
})

test_that("type validation rejects unknown types", {
  expect_error(nongenetic_feature("height", 1, db_or_freqs = c(a = 1)),
               "`type` must be one of")
  expect_error(nongenetic_feature(NA_character_, 1),
               "single non-empty character")
})

test_that("categorical db must be a normalised named vector", {
  expect_error(
    nongenetic_feature("sex", "F", db_or_freqs = c(F = 0.5, M = 0.4),
                       error = 0.05),
    "must sum to 1"
  )
  expect_error(
    nongenetic_feature("sex", "F", db_or_freqs = c(0.5, 0.5),
                       error = 0.05),
    "must be named"
  )
  expect_error(
    nongenetic_feature("sex", "F", db_or_freqs = c(F = 1),
                       error = 0.05),
    "at least two categories"
  )
})

test_that("categorical scalar error must lie in [0, 1)", {
  expect_error(
    nongenetic_feature("sex", "F", db_or_freqs = c(F = .5, M = .5),
                       error = 1),
    "finite scalar in \\[0, 1\\)"
  )
})

test_that("observed must match a declared category", {
  expect_error(
    nongenetic_feature("sex", "X", db_or_freqs = c(F = .5, M = .5),
                       error = 0.05),
    "is not one of the feature categories"
  )
})

test_that("region feature works with a migration confusion matrix", {
  cats <- c("N", "C", "S")
  E <- matrix(c(0.9, 0.05, 0.05,
                0.1, 0.8, 0.1,
                0.05, 0.15, 0.8), nrow = 3, byrow = TRUE,
              dimnames = list(cats, cats))
  f <- nongenetic_feature(
    type = "region",
    observed = "C",
    db_or_freqs = c(N = 0.5, C = 0.3, S = 0.2),
    error = E
  )
  expect_identical(f$feature_class, "categorical")
  expect_true(is.matrix(f$error))
  expect_identical(f$observed, "C")
})

test_that("print method runs for every feature class", {
  f_cat <- nongenetic_feature("sex", "F", db_or_freqs = c(F = .5, M = .5),
                              error = 0.05)
  f_con <- nongenetic_feature("age", 40, error = 0.05)
  f_dat <- nongenetic_feature("birthdate", 30,
                              error = c(1, 4, 60, 11, 6, 4, 4))
  expect_output(print(f_cat), "<nongenetic_feature>")
  expect_output(print(f_con), "continuous")
  expect_output(print(f_dat), "discrepancy bins")
  expect_identical(print(f_cat), f_cat)
})
