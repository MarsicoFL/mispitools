# compute_reference_prop() gives the population frequency of each
# hair / skin / eye combination. On a hand-made table we know the answer.

test_that("compute_reference_prop matches hand counts", {
  # 5 people: (1,1,1) twice, (2,1,1), (2,2,1), (2,2,2) once each
  d <- data.frame(hair_colour = c(1, 1, 2, 2, 2),
                  skin_colour = c(1, 1, 1, 2, 2),
                  eye_colour  = c(1, 1, 1, 1, 2))
  out <- compute_reference_prop(d)

  expect_named(out, c("hair_colour", "skin_colour", "eye_colour", "f_h_s_y"))
  # every combination of the observed levels is listed, 2 x 2 x 2
  expect_equal(nrow(out), 8)
  expect_equal(sum(out$f_h_s_y), 1)

  key <- paste(out$hair_colour, out$skin_colour, out$eye_colour)
  f <- setNames(out$f_h_s_y, key)
  expect_equal(f[["1 1 1"]], 2 / 5)
  expect_equal(f[["2 1 1"]], 1 / 5)
  expect_equal(f[["2 2 1"]], 1 / 5)
  expect_equal(f[["2 2 2"]], 1 / 5)
  # combinations nobody has come out as 0, not missing
  expect_equal(f[["1 2 2"]], 0)
})

test_that("compute_reference_prop agrees with sim_reference_pop counts", {
  pop <- sim_reference_pop(n = 300, seed = 1)
  out <- compute_reference_prop(pop)

  expect_equal(sum(out$f_h_s_y), 1)
  # spot check one combination against a direct count
  i <- which.max(out$f_h_s_y)
  n_direct <- sum(pop$hair_colour == out$hair_colour[i] &
                  pop$skin_colour == out$skin_colour[i] &
                  pop$eye_colour  == out$eye_colour[i])
  expect_equal(out$f_h_s_y[i], n_direct / 300)
})

test_that("compute_reference_prop rejects bad input", {
  expect_error(compute_reference_prop(list(a = 1)), "data.frame")
  expect_error(compute_reference_prop(data.frame(hair_colour = 1)),
               "skin_colour, eye_colour")
  empty <- data.frame(hair_colour = numeric(0), skin_colour = numeric(0),
                      eye_colour = numeric(0))
  expect_error(compute_reference_prop(empty), "no rows")
})
