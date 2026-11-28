# The LR calculators behind mispitools_app(), checked against values
# worked out by hand.

test_that("app_lr_sex matches the hand calculation", {
  # match: (1 - eps) / P(sex); mismatch: eps / P(sex)
  expect_equal(unname(app_lr_sex("F", "F", 0.05, 0.5)), 0.95 / 0.5)
  expect_equal(unname(app_lr_sex("F", "M", 0.05, 0.5)), 0.05 / 0.5)
  expect_equal(unname(app_lr_sex("M", "M", 0.01, 0.51)), 0.99 / 0.49)
  expect_equal(unname(app_lr_sex("M", "F", 0.01, 0.51)), 0.01 / 0.51)
  # no error at all: a mismatch excludes
  expect_equal(unname(app_lr_sex("F", "M", 0, 0.5)), 0)
})

test_that("app_lr_age matches the hand calculation", {
  # 40 +/- 5: P(in range | H2) = 10 / 80
  expect_equal(app_lr_age(40, 5, 42, 0.05), 0.95 / (10 / 80))
  expect_equal(app_lr_age(40, 5, 60, 0.05), 0.05 / (70 / 80))
  # interval ends count as in range
  expect_equal(app_lr_age(40, 5, 35, 0.05), app_lr_age(40, 5, 45, 0.05))
  expect_equal(app_lr_age(40, 5, 45, 0.05), 0.95 / (10 / 80))
})

test_that("app_lr_age clamps the interval to ages 1 to 80", {
  # 78 +/- 6 becomes [72, 80]: P(in range | H2) = 8 / 80
  expect_equal(app_lr_age(78, 6, 79, 0.05), 0.95 / (8 / 80))
  # 82 is inside 78 +/- 6 but above 80, so it counts as out of range
  expect_equal(app_lr_age(78, 6, 82, 0.05), 0.05 / (72 / 80))
  # 3 +/- 5 becomes [1, 8]
  expect_equal(app_lr_age(3, 5, 2, 0.05), 0.95 / (7 / 80))
})

test_that("app_lr_age handles zero-width and full-width intervals", {
  # zero width: P(in range | H2) is floored at 0.001
  expect_equal(app_lr_age(40, 0, 40, 0.05), 0.95 / 0.001)
  # wider than the population: clamped to [1, 80], so 79 / 80. The 0.999
  # cap in the function is never reached, because of the clamping.
  expect_equal(app_lr_age(40, 60, 10, 0.05), 0.95 / (79 / 80))
})

test_that("app_lr_hair matches the hand calculation", {
  freq <- c(Black = 0.08, Brown = 0.45, Blonde = 0.30, Red = 0.07, Gray = 0.10)
  # match: (1 - epc) / freq
  expect_equal(unname(app_lr_hair(2, 2, 0.1, freq)), 0.9 / 0.45)
  expect_equal(unname(app_lr_hair(4, 4, 0.1, freq)), 0.9 / 0.07)
  # mismatch: the error mass is spread evenly over the other 4 colours
  expect_equal(unname(app_lr_hair(2, 4, 0.1, freq)), (0.1 / 4) / 0.07)
  expect_equal(unname(app_lr_hair(1, 5, 0.4, freq)), (0.4 / 4) / 0.10)
})

test_that("app_lr_hair LRs average to 1 under H2", {
  # sum over observed colours of P(obs | H2) * LR = sum of P(obs | H1) = 1
  freq <- c(0.08, 0.45, 0.30, 0.07, 0.10)
  for (mpc in 1:5) {
    lr <- vapply(1:5, function(o) app_lr_hair(mpc, o, 0.1, freq), numeric(1))
    expect_equal(sum(freq * lr), 1)
  }
})
