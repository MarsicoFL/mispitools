test_that("migrants preliminary simulations preserve numeric columns", {
  db_poi <- sim_poi_prelim(casetype = "migrants", numsims = 100, seed = 456)
  expect_named(db_poi, c("UHR-ID", "Age", "Sex", "Height", "Region"))
  expect_type(db_poi[["UHR-ID"]], "double")
  expect_type(db_poi$Age, "integer")
  expect_type(db_poi$Height, "double")
  expect_type(db_poi$Sex, "character")
  expect_type(db_poi$Region, "character")

  db_mp <- sim_mp_prelim(casetype = "migrants", numsims = 100, seed = 456)
  expect_named(db_mp, c("UHR-ID", "Age", "Sex", "Height", "Region"))
  expect_type(db_mp[["UHR-ID"]], "double")
  expect_type(db_mp$Age, "integer")
  expect_type(db_mp$Height, "double")
  expect_type(db_mp$Sex, "character")
  expect_type(db_mp$Region, "character")
})

test_that("sim_lr_prelim uses numeric comparisons on migrants database", {
  db_poi <- sim_poi_prelim(casetype = "migrants", numsims = 100, seed = 456)
  n_match <- sum(db_poi$Age < 14 & db_poi$Age > 4)
  expect_gt(n_match, 0)

  expect_no_warning(
    res <- sim_lr_prelim(
      "age",
      numsims = 50,
      int = 5,
      ErrorRate = 0.05,
      MP = 9,
      database = db_poi
    )
  )
  expected_match_lr <- (1 - 0.05) / (n_match / nrow(db_poi))
  expected_mismatch_lr <- 0.05 / (1 - n_match / nrow(db_poi))
  expect_true(all(res$Related %in% c(expected_match_lr, expected_mismatch_lr)))
  expect_true(all(res$Unrelated %in% c(expected_match_lr, expected_mismatch_lr)))
})
