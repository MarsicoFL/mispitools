## F2.3 — mutation_matrix_cpp() vs mutation_matrix_R() and analytic.
##
## Validates that the pure-C++ builders for None / Equal / Stepwise
## reproduce the R reference engine bit-for-bit (1e-12 tol) and match
## closed-form analytic targets.

mat_cpp <- function(K, kind, rate = 0.0, range = 0.0,
                    numeric_labels = numeric(0), afreq = numeric(0)) {
  mispitools:::mutation_matrix_cpp(
    K = as.integer(K),
    mutation_kind = as.integer(kind),
    mutation_rate = as.double(rate),
    mutation_range = as.double(range),
    numeric_labels = as.numeric(numeric_labels),
    afreq = as.numeric(afreq)
  )
}

mat_R <- function(model, K, alleles = NULL) {
  mispitools:::mutation_matrix_R(model, K = K, alleles = alleles)
}

# ---------------------------------------------------------------------------
# None — identity
# ---------------------------------------------------------------------------

test_that("None: K=3 returns identity", {
  M <- mat_cpp(K = 3L, kind = 0L)
  expect_equal(M, diag(3), tolerance = 0)
})

test_that("None: matches mutation_matrix_R for K=2..5", {
  for (K in 2:5) {
    mod <- list(model = "none", rate = 0)
    expect_equal(mat_cpp(K, 0L), mat_R(mod, K), tolerance = 1e-15,
                 info = paste0("None K=", K))
  }
})

test_that("None: K=0 raises an error", {
  expect_error(mat_cpp(K = 0L, kind = 0L), "must be positive")
})

# ---------------------------------------------------------------------------
# Equal — analytic + R-ref cross-check
# ---------------------------------------------------------------------------

test_that("Equal K=2: closed-form (1-R) on diag, R off", {
  R <- 0.005
  M <- mat_cpp(K = 2L, kind = 1L, rate = R)
  expected <- matrix(c(1 - R, R, R, 1 - R), nrow = 2L, byrow = TRUE)
  expect_equal(M, expected, tolerance = 1e-15)
})

test_that("Equal K=4: off = R/(K-1) on every j != i", {
  R <- 0.01
  K <- 4L
  M <- mat_cpp(K, 1L, rate = R)
  off <- R / (K - 1)
  expected <- matrix(off, nrow = K, ncol = K)
  diag(expected) <- 1 - R
  expect_equal(M, expected, tolerance = 1e-15)
})

test_that("Equal: rows sum to 1 for K=2..8 and several rates", {
  for (K in c(2L, 3L, 5L, 8L)) {
    for (R in c(0, 1e-6, 0.001, 0.05, 0.5, 1)) {
      M <- mat_cpp(K, 1L, rate = R)
      expect_equal(rowSums(M), rep(1, K), tolerance = 1e-15,
                   info = paste0("Equal K=", K, " R=", R))
      expect_true(all(M >= -1e-15 & M <= 1 + 1e-15),
                  info = paste0("Equal K=", K, " R=", R))
    }
  }
})

test_that("Equal: matches mutation_matrix_R bit-for-bit", {
  for (K in 2:6) {
    for (R in c(0, 1e-4, 0.005, 0.1, 0.999)) {
      mod <- list(model = "equal", rate = R)
      expect_equal(mat_cpp(K, 1L, rate = R), mat_R(mod, K),
                   tolerance = 1e-15,
                   info = paste0("Equal K=", K, " R=", R))
    }
  }
})

test_that("Equal: rate=0 collapses to identity", {
  M <- mat_cpp(K = 5L, kind = 1L, rate = 0)
  expect_equal(M, diag(5), tolerance = 1e-15)
})

test_that("Equal: validates K and rate", {
  expect_error(mat_cpp(K = 1L, kind = 1L, rate = 0.005), "K >= 2")
  expect_error(mat_cpp(K = 3L, kind = 1L, rate = -0.001), "\\[0, 1\\]")
  expect_error(mat_cpp(K = 3L, kind = 1L, rate = 1.0001), "\\[0, 1\\]")
})

# ---------------------------------------------------------------------------
# Stepwise — analytic + R-ref cross-check + boundary cases
# ---------------------------------------------------------------------------

test_that("Stepwise K=2: closed form independent of range", {
  ## With two alleles, |s_j - s_i| > 0 for the single off-diagonal entry,
  ## so the normalized weight is 1 regardless of `range`. Result is the
  ## same as Equal-rate for K = 2.
  labs <- c(10, 11)
  R <- 0.003
  for (r in c(0.05, 0.1, 0.5)) {
    M <- mat_cpp(K = 2L, kind = 2L, rate = R, range = r,
                 numeric_labels = labs)
    expected <- matrix(c(1 - R, R, R, 1 - R), nrow = 2L, byrow = TRUE)
    expect_equal(M, expected, tolerance = 1e-15,
                 info = paste0("Stepwise K=2 r=", r))
  }
})

test_that("Stepwise K=3 unit-spaced: analytic closed form", {
  R <- 0.004
  r <- 0.1
  labs <- c(10, 11, 12)
  M <- mat_cpp(K = 3L, kind = 2L, rate = R, range = r,
               numeric_labels = labs)

  ## Row 0 (parent allele 10): distances {-, 1, 2} -> w = (0, r, r^2)
  ##   sw_0 = r + r^2; off = R / sw_0 * w_j
  sw0 <- r + r^2
  expect_equal(M[1, 1], 1 - R, tolerance = 1e-15)
  expect_equal(M[1, 2], (R / sw0) * r, tolerance = 1e-15)
  expect_equal(M[1, 3], (R / sw0) * r^2, tolerance = 1e-15)
  ## Row 1 (parent 11): distances {1, -, 1} -> w = (r, 0, r); sw_1 = 2r
  expect_equal(M[2, 1], R / 2, tolerance = 1e-15)
  expect_equal(M[2, 2], 1 - R, tolerance = 1e-15)
  expect_equal(M[2, 3], R / 2, tolerance = 1e-15)
  ## Row 2 mirrors row 0
  expect_equal(M[3, 1], (R / sw0) * r^2, tolerance = 1e-15)
  expect_equal(M[3, 2], (R / sw0) * r, tolerance = 1e-15)
  expect_equal(M[3, 3], 1 - R, tolerance = 1e-15)
  expect_equal(rowSums(M), rep(1, 3), tolerance = 1e-15)
})

test_that("Stepwise: matches mutation_matrix_R bit-for-bit", {
  cases <- list(
    list(alleles = c("10", "11", "12", "13"), R = 0.005, r = 0.1),
    list(alleles = c("8", "9.3", "10", "11.3"), R = 0.001, r = 0.5),
    list(alleles = c("12", "14", "16", "18", "20"), R = 0.02, r = 0.25),
    list(alleles = c("15", "16"), R = 0.5, r = 0.1)
  )
  for (cs in cases) {
    K <- length(cs$alleles)
    mod <- list(model = "stepwise", rate = cs$R, ratio = cs$r)
    labs <- as.numeric(cs$alleles)
    M_cpp <- mat_cpp(K, 2L, rate = cs$R, range = cs$r, numeric_labels = labs)
    M_R <- mat_R(mod, K = K, alleles = cs$alleles)
    expect_equal(M_cpp, M_R, tolerance = 1e-15,
                 info = paste0("Stepwise alleles=",
                               paste(cs$alleles, collapse = ","),
                               " R=", cs$R, " r=", cs$r))
  }
})

test_that("Stepwise: rows sum to 1 and entries are in [0, 1]", {
  R <- 0.01
  r <- 0.1
  labs <- c(10, 11, 12, 13, 14, 15)
  M <- mat_cpp(K = 6L, kind = 2L, rate = R, range = r, numeric_labels = labs)
  expect_equal(rowSums(M), rep(1, 6), tolerance = 1e-15)
  expect_true(all(M >= 0 & M <= 1))
})

test_that("Stepwise: continuity R->0 collapses to identity", {
  labs <- c(10, 11, 12)
  M <- mat_cpp(K = 3L, kind = 2L, rate = 0, range = 0.1,
               numeric_labels = labs)
  expect_equal(M, diag(3), tolerance = 1e-15)
})

test_that("Stepwise: range = 0 raises (all off-diag weights zero)", {
  labs <- c(10, 11, 12)
  expect_error(
    mat_cpp(K = 3L, kind = 2L, rate = 0.005, range = 0,
            numeric_labels = labs),
    "sum to zero"
  )
})

test_that("Stepwise: validates K, rate, range, labels", {
  expect_error(
    mat_cpp(K = 1L, kind = 2L, rate = 0.005, range = 0.1,
            numeric_labels = c(10)),
    "K >= 2"
  )
  expect_error(
    mat_cpp(K = 3L, kind = 2L, rate = -0.01, range = 0.1,
            numeric_labels = c(10, 11, 12)),
    "\\[0, 1\\]"
  )
  expect_error(
    mat_cpp(K = 3L, kind = 2L, rate = 0.005, range = -0.1,
            numeric_labels = c(10, 11, 12)),
    "range must be finite and >= 0"
  )
  ## numeric_labels.size() != K
  expect_error(
    mat_cpp(K = 3L, kind = 2L, rate = 0.005, range = 0.1,
            numeric_labels = c(10, 11)),
    "numeric_labels.size\\(\\) must equal K"
  )
  ## non-finite label
  expect_error(
    mat_cpp(K = 3L, kind = 2L, rate = 0.005, range = 0.1,
            numeric_labels = c(10, NaN, 12)),
    "every numeric_label must be"
  )
})

# ---------------------------------------------------------------------------
# Dispatcher edges
# ---------------------------------------------------------------------------

test_that("Dispatcher: Proportional still returns an informative error", {
  expect_error(
    mat_cpp(K = 3L, kind = 3L, rate = 0.005),
    "Proportional model is not implemented"
  )
})

test_that("Dispatcher: Asymmetric without afreq is rejected", {
  expect_error(
    mat_cpp(K = 3L, kind = 4L, rate = 0.005, range = 0.1),
    "Asymmetric \\(Dawid\\) model requires allele frequencies"
  )
})

test_that("Dispatcher: out-of-range kind is rejected at the binding", {
  expect_error(
    mat_cpp(K = 3L, kind = 5L),
    "mutation_kind out of range"
  )
})
