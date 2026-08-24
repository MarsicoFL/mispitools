## F0.5 — R -> Rcpp -> mispitools::core::*_placeholder roundtrip.
##
## Verifies that:
##   * the package shared object links and loads,
##   * the .Call() registration table in src/RcppExports.cpp matches the
##     stubs in R/RcppExports.R,
##   * integer marshalling works in both directions across the Rcpp
##     boundary, and
##   * every core/ placeholder returns x + 1 as F0.4 fixed.
##
## Each placeholder is exercised once with x = 41 (canonical "answer to
## life" smoke test from the F0.4 standalone driver) and once with a
## sweep of small/negative/large inputs.

placeholder_bindings <- list(
  cpp_pedigree_placeholder,
  cpp_marker_placeholder,
  cpp_mutation_models_placeholder,
  cpp_cpt_engine_placeholder,
  cpp_kl_engine_placeholder,
  cpp_lr_dist_placeholder,
  cpp_nongenetic_lr_placeholder,
  cpp_evidence_combine_placeholder,
  cpp_concentration_placeholder,
  cpp_decision_placeholder,
  cpp_linkage_placeholder
)

test_that("all 11 core placeholders return 42 when called with 41", {
  results <- vapply(placeholder_bindings, function(fn) fn(41L), integer(1))
  expect_equal(results, rep(42L, length(placeholder_bindings)))
})

test_that("placeholders behave as add_one over a range of integer inputs", {
  inputs <- c(-100L, -1L, 0L, 1L, 7L, 99L, 10000L)
  for (fn in placeholder_bindings) {
    expect_equal(vapply(inputs, fn, integer(1)), inputs + 1L)
  }
})

test_that("each binding is registered with one argument", {
  syms <- c(
    "_mispitools_cpp_pedigree_placeholder",
    "_mispitools_cpp_marker_placeholder",
    "_mispitools_cpp_mutation_models_placeholder",
    "_mispitools_cpp_cpt_engine_placeholder",
    "_mispitools_cpp_kl_engine_placeholder",
    "_mispitools_cpp_lr_dist_placeholder",
    "_mispitools_cpp_nongenetic_lr_placeholder",
    "_mispitools_cpp_evidence_combine_placeholder",
    "_mispitools_cpp_concentration_placeholder",
    "_mispitools_cpp_decision_placeholder",
    "_mispitools_cpp_linkage_placeholder"
  )
  routines <- getDLLRegisteredRoutines("mispitools", addNames = TRUE)
  call_routines <- routines[[".Call"]]
  expect_true(all(syms %in% names(call_routines)))
  for (s in syms) {
    expect_equal(call_routines[[s]]$numParameters, 1L,
                 info = paste("routine", s, "should take 1 SEXP argument"))
  }
})
