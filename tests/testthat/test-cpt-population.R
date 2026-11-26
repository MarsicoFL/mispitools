# cpt_population() is P(sex, age group, hair colour | H2) as a product of
# independent marginals, so every cell can be written down by hand.
# It is soft-deprecated, so the one-time notice is silenced here.

cpt_pop <- function(...) suppressMessages(cpt_population(...))

test_that("cpt_population default table matches the hand calculation", {
  out <- cpt_pop()
  # MPa = 40, MPr = 6: P(T1) = 12 / 80 = 0.15
  sex_age <- c(`F-T1` = 0.5 * 0.15, `F-T0` = 0.5 * 0.85,
               `M-T1` = 0.5 * 0.15, `M-T0` = 0.5 * 0.85)
  hair <- c(0.3, 0.2, 0.25, 0.15, 0.1)
  expected <- outer(sex_age, hair)
  colnames(expected) <- as.character(1:5)

  expect_equal(out, expected)
  expect_equal(out["F-T1", "1"], 0.0225)
  expect_equal(out["M-T0", "5"], 0.0425)
})

test_that("cpt_population margins recover the inputs", {
  propS <- c(0.6, 0.4)
  propC <- c(0.4, 0.25, 0.2, 0.1, 0.05)
  out <- cpt_pop(propS = propS, MPa = 35, MPr = 3, propC = propC)
  t1 <- 6 / 80

  expect_equal(sum(out), 1)
  expect_equal(unname(colSums(out)), propC)
  expect_equal(unname(rowSums(out)),
               c(0.6 * t1, 0.6 * (1 - t1), 0.4 * t1, 0.4 * (1 - t1)))
})

test_that("cpt_population normalises proportions that don't sum to 1", {
  expect_warning(a <- cpt_pop(propS = c(1, 1)), "propS does not sum to 1")
  expect_equal(a, cpt_pop())
  expect_warning(b <- cpt_pop(propC = c(3, 2, 2.5, 1.5, 1)),
                 "propC does not sum to 1")
  expect_equal(b, cpt_pop())
})

test_that("cpt_population rejects bad input", {
  expect_error(cpt_pop(propS = 1), "length 2")
  expect_error(cpt_pop(propS = c(-0.5, 1.5)), "between 0 and 1")
  expect_error(cpt_pop(MPa = 0), "between 1 and 80")
  expect_error(cpt_pop(MPr = -1), "non-negative")
  expect_error(cpt_pop(propC = c(0.5, 0.5)), "length 5")
  expect_error(cpt_pop(propC = c(-0.1, 0.3, 0.3, 0.3, 0.2)), "non-negative")
})
