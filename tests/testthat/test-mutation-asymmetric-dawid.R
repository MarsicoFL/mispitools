## F5.1 — core::mutation_matrix_asymmetric() = Dawid (2002) unified
## asymmetric / reversible model.
##
## Oracle: pedmut::mutationMatrix(model = "dawid", ...) at 1e-8
## (entry-wise), per SCOUT_pedmut_asymmetric.md. Plus a deterministic
## closed-form anchor, reversibility / row-stochastic invariants, the
## undefined-model rate cap, and input validation.
##
## Note on the misnomer: STATE.md / ROADMAP label this hito "Familias
## bias u". SCOUT F4.6 established there is no directional bias `u` in
## Familias or pedmut (the candidate src/alsys.cpp block is dead,
## commented-out code). The canonical asymmetric model of the pedsuite
## ecosystem is Dawid 2002; that is what is implemented here. A true
## directional-slippage bias has no oracle and is intentionally NOT
## implemented in 2.0.

dawid_cpp <- function(afreq, rate, range) {
  mispitools:::mutation_matrix_cpp(
    K = length(afreq),
    mutation_kind = 4L,
    mutation_rate = as.double(rate),
    mutation_range = as.double(range),
    afreq = as.numeric(afreq)
  )
}

# ---------------------------------------------------------------------------
# Deterministic closed-form anchor (independent of pedmut). Hand-verified
# against pedmut::mutationMatrix("dawid") for afreq = c(.2,.5,.3),
# rate = 0.002, range = 0.1: see SCOUT_pedmut_asymmetric.md algorithm.
# ---------------------------------------------------------------------------

test_that("Dawid K=3 reproduces the hand-computed matrix", {
  M <- dawid_cpp(c(0.2, 0.5, 0.3), rate = 0.002, range = 0.1)
  expected <- rbind(
    c(0.9973809523809524, 0.0023809523809524, 0.0002380952380952),
    c(0.0009523809523810, 0.9980952380952381, 0.0009523809523810),
    c(0.0001587301587302, 0.0015873015873016, 0.9982539682539683)
  )
  expect_equal(M, expected, tolerance = 1e-12)
})

test_that("Dawid rows are stochastic and the matrix is reversible", {
  p <- c(0.2, 0.5, 0.3)
  M <- dawid_cpp(p, rate = 0.003, range = 0.2)
  expect_equal(rowSums(M), rep(1, 3), tolerance = 1e-14)
  # Detailed balance: p[i] M[i,j] == p[j] M[j,i].
  db <- outer(p, rep(1, 3)) * M
  expect_equal(db, t(db), tolerance = 1e-14)
})

# ---------------------------------------------------------------------------
# Oracle: pedmut::mutationMatrix(model = "dawid"). Five Argentina
# integer-allele STR markers x {rate} x range = 0.1, tol 1e-8.
# ---------------------------------------------------------------------------

DAWID_MARKERS <- c("D8S1179", "D7S820", "CSF1PO", "D3S1358", "D13S317")

for (mk in DAWID_MARKERS) {
  for (rt in c(0.001, 0.005)) {
    local({
      mk_ <- mk; rt_ <- rt
      test_that(sprintf("Dawid vs pedmut: %s rate=%g", mk_, rt_), {
        skip_if_not_installed("pedmut")
        data(Argentina, package = "mispitools", envir = environment())
        f <- top_k_freqs(Argentina, mk_, top_k = 5L)
        M_cpp <- dawid_cpp(as.numeric(f), rate = rt_, range = 0.1)
        M_ref <- unclass(pedmut::mutationMatrix(
          model = "dawid", alleles = names(f),
          afreq = as.numeric(f), rate = rt_, range = 0.1))
        attributes(M_ref) <- list(dim = dim(M_ref))
        expect_equal(M_cpp, M_ref, tolerance = 1e-8,
                     info = sprintf("%s rate=%g", mk_, rt_))
        expect_true(pedmut::isReversible(
          pedmut::mutationMatrix("dawid", alleles = names(f),
                                 afreq = as.numeric(f),
                                 rate = rt_, range = 0.1)))
      })
    })
  }
}

# ---------------------------------------------------------------------------
# Undefined-model rate cap (pedmut::maxRate UW bound).
# ---------------------------------------------------------------------------

test_that("dawid_max_rate_cpp matches the analytic UW cap", {
  cap <- mispitools:::dawid_max_rate_cpp(c(0.2, 0.5, 0.3), range = 0.1)
  expect_equal(cap, 0.763636363636364, tolerance = 1e-10)
})

test_that("rate above the cap yields an undefined-model error with the cap", {
  expect_error(
    dawid_cpp(c(0.2, 0.5, 0.3), rate = 0.9, range = 0.1),
    "Dawid model undefined; max rate for the given input is 0.7636"
  )
})

test_that("cap agrees with pedmut's reported maximum rate", {
  skip_if_not_installed("pedmut")
  p <- c(0.2, 0.5, 0.3)
  cap <- mispitools:::dawid_max_rate_cpp(p, range = 0.1)
  msg <- tryCatch(
    pedmut::mutationMatrix("dawid", alleles = c("12", "14", "15"),
                           afreq = p, rate = 0.9, range = 0.1),
    error = function(e) conditionMessage(e))
  ref <- as.numeric(sub(".*is: *", "", msg))
  expect_equal(cap, ref, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Input validation.
# ---------------------------------------------------------------------------

test_that("Dawid rejects K < 2", {
  expect_error(dawid_cpp(c(1.0), rate = 0.001, range = 0.1),
               "requires K >= 2")
})

test_that("Dawid rejects range outside (0, 1)", {
  expect_error(dawid_cpp(c(0.5, 0.5), rate = 0.001, range = 0),
               "range finite in \\(0, 1\\)")
  expect_error(dawid_cpp(c(0.5, 0.5), rate = 0.001, range = 1),
               "range finite in \\(0, 1\\)")
})

test_that("Dawid rejects non-positive or non-finite frequencies", {
  expect_error(dawid_cpp(c(0.5, 0.0, 0.5), rate = 0.001, range = 0.1),
               "finite and strictly positive")
})

test_that("Dawid rejects an out-of-range rate", {
  expect_error(dawid_cpp(c(0.5, 0.5), rate = 1.5, range = 0.1),
               "rate must be finite in \\[0, 1\\]")
})

test_that("dawid_max_rate_cpp validates its inputs", {
  expect_error(mispitools:::dawid_max_rate_cpp(numeric(0), range = 0.1),
               "requires K >= 2")
  expect_error(mispitools:::dawid_max_rate_cpp(c(0.5, 0.5), range = 1),
               "range finite in \\(0, 1\\)")
})
