## F3.4 — N-marker batch + mutation matrix cache.
##
## Two layers under test:
##   * cpp_per_marker_kl_batch() — the C++ binding called directly with a
##     flattened pedigree + per-marker freqs/labels/mutation.
##   * per_marker_kl_profile() — the public R API, which detects shared
##     pedigree topology and routes to the batch path; otherwise falls
##     back to the scalar loop landed in F3.2.
##
## All numerical checks are bit-exact (1e-15) versus the scalar entry
## points landed in F3.1 / F3.2; the batch path must not perturb the
## per-marker outputs in any way. The cache_hits diagnostic is verified
## directly to confirm the matrix is reused across markers with matching
## (kind, K, rate, range, labels) signatures.

skip_if_no_pedtools <- function() {
  testthat::skip_if_not_installed("pedtools")
}

batch_cols <- c("marker", "e_log10_lr_h1", "e_log10_lr_h2",
                "kl_h1_to_h2", "kl_h2_to_h1",
                "abs_cont_violations_h1", "abs_cont_violations_h2",
                "mass_violations_h1", "mass_violations_h2")

# ---------------------------------------------------------------------------
# Batch matches the scalar loop bit-for-bit
# ---------------------------------------------------------------------------

test_that("per_marker_kl_profile batch == scalar loop (mut=equal, K=2/K=3)", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  models <- list(
    M1 = marker_model(ped, "M1", c("a" = 0.4, "b" = 0.6),
                      mutation = list(model = "equal", rate = 0.005)),
    M2 = marker_model(ped, "M2", c("a" = 0.2, "b" = 0.3, "c" = 0.5),
                      mutation = list(model = "equal", rate = 0.005))
  )
  out <- per_marker_kl_profile(models)
  expect_named(out, batch_cols)
  expect_equal(nrow(out), 2L)
  for (i in seq_along(models)) {
    s <- per_marker_kl(models[[i]])
    for (col in setdiff(batch_cols, "marker")) {
      expect_equal(out[[col]][i], s[[col]],
                   tolerance = 1e-15, info = paste(col, i))
    }
  }
})

test_that("per_marker_kl_profile batch == scalar loop (stepwise + halfSibPed)", {
  skip_if_no_pedtools()
  ped <- pedtools::halfSibPed()
  models <- list(
    marker_model(ped, "M1",
                 c("12" = 0.3, "13" = 0.4, "14" = 0.3),
                 mutation = list(model = "stepwise", rate = 0.005, ratio = 0.1)),
    marker_model(ped, "M2",
                 c("15" = 0.2, "16" = 0.3, "17" = 0.5),
                 mutation = list(model = "stepwise", rate = 0.005, ratio = 0.1)),
    marker_model(ped, "M3",
                 c("a" = 0.4, "b" = 0.6),
                 mutation = list(model = "equal", rate = 0.01))
  )
  out <- per_marker_kl_profile(models)
  for (i in seq_along(models)) {
    s <- per_marker_kl(models[[i]])
    for (col in setdiff(batch_cols, "marker")) {
      expect_equal(out[[col]][i], s[[col]],
                   tolerance = 1e-15, info = paste(col, i))
    }
  }
})

test_that("per_marker_kl_profile batch preserves +Inf branch (nuclearPed(2), mut=none)", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(2)
  freqs <- c("a" = 0.4, "b" = 0.6)
  models <- list(
    marker_model(ped, "MA", freqs),
    marker_model(ped, "MB", freqs)
  )
  out <- per_marker_kl_profile(models)
  for (i in seq_along(models)) {
    s <- per_marker_kl(models[[i]])
    expect_identical(out$kl_h2_to_h1[i], s$kl_h2_to_h1, info = i)
    expect_identical(out$e_log10_lr_h2[i], s$e_log10_lr_h2, info = i)
    expect_equal(out$kl_h1_to_h2[i], s$kl_h1_to_h2,
                 tolerance = 1e-15, info = i)
    expect_gte(out$abs_cont_violations_h1[i], 1L)
    expect_gt(out$mass_violations_h1[i], 0)
  }
})

# ---------------------------------------------------------------------------
# Mutation matrix cache diagnostics
# ---------------------------------------------------------------------------

test_that("batch caches the mutation matrix when (kind, K, rate, range) matches", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  ## Three K=2 + 1 K=3 markers, all with mut=equal rate=0.005. Expect
  ## two distinct cache entries (K=2 and K=3) → 2 misses, 2 hits.
  models <- list(
    marker_model(ped, "A1", c("a" = 0.4, "b" = 0.6),
                 mutation = list(model = "equal", rate = 0.005)),
    marker_model(ped, "A2", c("p" = 0.7, "q" = 0.3),
                 mutation = list(model = "equal", rate = 0.005)),
    marker_model(ped, "B1", c("a" = 0.2, "b" = 0.3, "c" = 0.5),
                 mutation = list(model = "equal", rate = 0.005)),
    marker_model(ped, "A3", c("a" = 0.5, "b" = 0.5),
                 mutation = list(model = "equal", rate = 0.005))
  )
  out <- per_marker_kl_profile(models)
  expect_identical(attr(out, "cache_hits"),   2L)
  expect_identical(attr(out, "cache_misses"), 2L)
})

test_that("batch caches mutation matrices for stepwise across identical labels", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  models <- list(
    marker_model(ped, "S1", c("12" = 0.3, "13" = 0.4, "14" = 0.3),
                 mutation = list(model = "stepwise", rate = 0.005, ratio = 0.1)),
    marker_model(ped, "S2", c("12" = 0.5, "13" = 0.2, "14" = 0.3),
                 mutation = list(model = "stepwise", rate = 0.005, ratio = 0.1))
  )
  out <- per_marker_kl_profile(models)
  expect_identical(attr(out, "cache_hits"),   1L)
  expect_identical(attr(out, "cache_misses"), 1L)
})

test_that("batch does NOT share cache across mismatched (kind, K, rate, range)", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  models <- list(
    marker_model(ped, "A", c("a" = 0.4, "b" = 0.6),
                 mutation = list(model = "equal", rate = 0.005)),
    marker_model(ped, "B", c("a" = 0.4, "b" = 0.6),
                 mutation = list(model = "equal", rate = 0.01)),
    marker_model(ped, "C", c("12" = 0.5, "13" = 0.5),
                 mutation = list(model = "stepwise", rate = 0.005,
                                 ratio = 0.1))
  )
  out <- per_marker_kl_profile(models)
  expect_identical(attr(out, "cache_hits"),   0L)
  expect_identical(attr(out, "cache_misses"), 3L)
})

# ---------------------------------------------------------------------------
# F3.4b: cached batch == uncached scalar route, bit-for-bit (1e-12)
# ---------------------------------------------------------------------------
#
# The scalar `per_marker_kl(model)` rebuilds the K x K mutation matrix
# independently for every marker (no cache), so it is the uncached oracle
# for the cached batch path. A profile with a repeated mutation signature
# exercises a real cache hit (asserted) while the numbers must remain
# identical to the per-marker rebuild down to 1e-12.

test_that("cached batch result is bit-for-bit the uncached scalar route", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  ## Markers 1/2/4 share (equal, K=2, rate=0.005) → 2 cache hits; marker 3
  ## is K=3 (fresh miss); marker 5 is stepwise (fresh miss). Net: >= 1 hit.
  models <- list(
    marker_model(ped, "G1", c("a" = 0.4, "b" = 0.6),
                 mutation = list(model = "equal", rate = 0.005)),
    marker_model(ped, "G2", c("p" = 0.8, "q" = 0.2),
                 mutation = list(model = "equal", rate = 0.005)),
    marker_model(ped, "G3", c("a" = 0.2, "b" = 0.3, "c" = 0.5),
                 mutation = list(model = "equal", rate = 0.005)),
    marker_model(ped, "G4", c("a" = 0.55, "b" = 0.45),
                 mutation = list(model = "equal", rate = 0.005)),
    marker_model(ped, "G5", c("12" = 0.3, "13" = 0.4, "14" = 0.3),
                 mutation = list(model = "stepwise", rate = 0.005,
                                 ratio = 0.1))
  )
  out <- per_marker_kl_profile(models)
  expect_gte(attr(out, "cache_hits"), 1L)
  for (i in seq_along(models)) {
    uncached <- per_marker_kl(models[[i]])
    for (col in setdiff(batch_cols, "marker")) {
      expect_equal(out[[col]][i], uncached[[col]],
                   tolerance = 1e-12, info = paste(col, i))
    }
  }
})

# ---------------------------------------------------------------------------
# Fallback path: heterogeneous topology must produce the same numbers
# ---------------------------------------------------------------------------

test_that("profile falls back to scalar loop when pedigrees differ", {
  skip_if_no_pedtools()
  models <- list(
    marker_model(pedtools::nuclearPed(1), "M1",
                 c("a" = 0.4, "b" = 0.6),
                 mutation = list(model = "equal", rate = 0.005)),
    marker_model(pedtools::halfSibPed(), "M2",
                 c("a" = 0.4, "b" = 0.6),
                 mutation = list(model = "equal", rate = 0.005))
  )
  out <- per_marker_kl_profile(models)
  expect_equal(nrow(out), 2L)
  expect_null(attr(out, "cache_hits"))
  for (i in seq_along(models)) {
    s <- per_marker_kl(models[[i]])
    for (col in setdiff(batch_cols, "marker")) {
      expect_equal(out[[col]][i], s[[col]],
                   tolerance = 1e-15, info = paste(col, i))
    }
  }
})

# ---------------------------------------------------------------------------
# Direct binding: shape + cache stats
# ---------------------------------------------------------------------------

test_that("cpp_per_marker_kl_batch returns the documented list shape", {
  skip_if_no_pedtools()
  ## Two K=2 markers with identical mutation parameters → matrix cached
  ## (1 miss + 1 hit). The third K=3 marker would be a fresh miss.
  freqs1 <- c(0.4, 0.6)
  freqs2 <- c(0.7, 0.3)
  res <- mispitools:::cpp_per_marker_kl_batch(
    father = c(-1L, -1L, 0L),
    mother = c(-1L, -1L, 1L),
    poi = 2L,
    freqs_list = list(freqs1, freqs2),
    mutation_kind = c(1L, 1L),
    mutation_rate = c(0.005, 0.005),
    mutation_range = c(0.0, 0.0),
    numeric_labels_list = list(rep(NA_real_, 2), rep(NA_real_, 2))
  )
  expect_named(res, c("e_log10_lr_h1", "e_log10_lr_h2",
                      "kl_h1_to_h2", "kl_h2_to_h1",
                      "abs_cont_violations_h1", "abs_cont_violations_h2",
                      "mass_violations_h1", "mass_violations_h2",
                      "cache_hits", "cache_misses"))
  expect_length(res$e_log10_lr_h1, 2L)
  expect_identical(res$cache_hits, 1L)
  expect_identical(res$cache_misses, 1L)
})

test_that("cpp_per_marker_kl_batch validates list lengths and kind range", {
  expect_error(
    mispitools:::cpp_per_marker_kl_batch(
      father = c(-1L, -1L, 0L),
      mother = c(-1L, -1L, 1L),
      poi = 2L,
      freqs_list = list(c(0.5, 0.5)),
      mutation_kind = c(1L, 1L),
      mutation_rate = c(0.005, 0.005),
      mutation_range = c(0.0, 0.0),
      numeric_labels_list = list(rep(NA_real_, 2))
    ),
    "per-marker vectors"
  )
  expect_error(
    mispitools:::cpp_per_marker_kl_batch(
      father = c(-1L, -1L, 0L),
      mother = c(-1L, -1L, 1L),
      poi = 2L,
      freqs_list = list(c(0.5, 0.5)),
      mutation_kind = c(7L),
      mutation_rate = c(0.005),
      mutation_range = c(0.0),
      numeric_labels_list = list(rep(NA_real_, 2))
    ),
    "mutation_kind out of range"
  )
})
