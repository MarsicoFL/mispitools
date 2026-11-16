# threshold_rates() checked against the sex LR, where we can work out the
# answer by hand.
#
# With no MP data (open search) the sex LR is either 1.9 (sex matches) or
# 0.1 (it doesn't), using the default error rate of 0.05. A true match has
# the right sex 95% of the time, a random person 50% of the time. So with
# any threshold between 0.1 and 1.9 we should get FNR = 0.05 and FPR = 0.5.

test_that("threshold_rates gives exact rates on a known sample", {
  # 95 true matches at 1.9 and 5 at 0.1; unrelated split 50/50
  x <- data.frame(Related   = c(rep(1.9, 95), rep(0.1, 5)),
                  Unrelated = c(rep(1.9, 50), rep(0.1, 50)))
  r <- suppressMessages(threshold_rates(x, threshold = 1))

  expect_equal(r$FNR, 0.05)
  expect_equal(r$TPR, 0.95)
  expect_equal(r$FPR, 0.50)
  expect_equal(r$TNR, 0.50)
})

test_that("threshold_rates handles thresholds outside the LR range", {
  x <- data.frame(Related   = c(rep(1.9, 95), rep(0.1, 5)),
                  Unrelated = c(rep(1.9, 50), rep(0.1, 50)))

  # below every LR: everyone counts as a match
  below <- suppressMessages(threshold_rates(x, threshold = 0.05))
  expect_equal(c(below$FNR, below$FPR), c(0, 1))

  # above every LR: nobody does
  above <- suppressMessages(threshold_rates(x, threshold = 3))
  expect_equal(c(above$FNR, above$FPR), c(1, 0))
})

test_that("threshold_rates agrees with the hand calculation on simulated LRs", {
  for (e in c(0.01, 0.05, 0.2)) {
    x <- sim_lr_prelim("sex", numsims = 20000, ErrorRate = e, seed = 1)
    r <- suppressMessages(threshold_rates(x, threshold = 1))

    # simulated, so allow some noise
    expect_lt(abs(r$FNR - e), 0.01)
    expect_lt(abs(r$FPR - 0.5), 0.02)
  }
})
