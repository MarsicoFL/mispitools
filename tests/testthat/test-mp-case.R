
skip_if_no_pedtools <- function() {
  testthat::skip_if_not_installed("pedtools")
}

sex_feature <- function(observed = "F", eps = 0.05) {
  nongenetic_feature(
    type = "sex", observed = observed,
    db_or_freqs = c(F = 0.5, M = 0.5), error = eps
  )
}

test_that("mp_case builds genetic models from a wide frequency table", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  case <- mp_case(
    genetic = list(ped = ped, freqs = Argentina,
                   markers = c("D3S1358", "THO1"),
                   mutation = list(model = "equal", rate = 1e-3))
  )
  expect_s3_class(case, "mp_case")
  expect_equal(case$n_markers, 2L)
  expect_equal(case$n_features, 0L)
  expect_length(case$genetic$models, 2L)
  expect_true(all(vapply(case$genetic$models, inherits, logical(1),
                         "marker_model")))
  expect_equal(names(case$genetic$models), c("D3S1358", "THO1"))
  for (m in case$genetic$models) {
    expect_equal(sum(m$freqs), 1, tolerance = 1e-12)
    expect_true(all(m$freqs > 0))
  }
})

test_that("mp_case accepts non-genetic only and a bare feature", {
  case <- mp_case(non_genetic = sex_feature())
  expect_s3_class(case, "mp_case")
  expect_equal(case$n_markers, 0L)
  expect_equal(case$n_features, 1L)
  expect_null(case$genetic)
})

test_that("mp_case rejects empty evidence and bad inputs", {
  expect_error(mp_case(), "at least one source of evidence")
  expect_error(
    mp_case(non_genetic = list(1, 2)),
    "must be 'nongenetic_feature'")
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  expect_error(
    mp_case(genetic = list(ped = ped)),
    "freqs.*is required")
  expect_error(
    mp_case(genetic = list(freqs = Argentina)),
    "ped.*is required")
  expect_error(
    mp_case(genetic = list(ped = ped, freqs = Argentina,
                           markers = "NOPE")),
    "absent from")
})

test_that("markov_se requires a transition list of the right length", {
  expect_error(
    mp_case(non_genetic = list(sex_feature(), sex_feature("M")),
            combine = "markov_se"),
    "requires `transition`")
})

test_that("evaluate_evidence returns the documented contract", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  case <- mp_case(
    genetic = list(ped = ped, freqs = Argentina,
                   markers = c("D3S1358", "THO1"),
                   mutation = list(model = "equal", rate = 1e-3)),
    non_genetic = list(sex = sex_feature())
  )
  res <- evaluate_evidence(case)

  expect_s3_class(res, "evidence_result")
  expect_s3_class(res$per_feature_kl, "data.frame")
  expect_named(res$per_feature_kl,
               c("feature", "kind", "e_log10_lr_h1", "e_log10_lr_h2",
                 "kl_h1_to_h2", "kl_h2_to_h1"))
  expect_equal(nrow(res$per_feature_kl), 3L)
  expect_equal(res$per_feature_kl$kind,
               c("genetic", "genetic", "non_genetic"))
  expect_s3_class(res$combined_lr_dist, "lr_dist")
  expect_length(res$per_feature_lr_dist, 3L)
  expect_named(res$per_feature_lr_dist,
               c("D3S1358", "THO1", "sex"))
  expect_true(is.finite(res$expected_log10_lr_h1))
  expect_true(all(c("threshold", "fpr", "fnr", "auc") %in%
                  names(res$decision_quantities)))
  expect_true(res$fragility$concentration >= 0 &&
              res$fragility$concentration <= 1)
  expect_s3_class(res$fragility$leave_one_out, "data.frame")
})

test_that("independent combine == lr_distribution over genetic markers", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  spec_markers <- c("D3S1358", "THO1")
  case <- mp_case(
    genetic = list(ped = ped, freqs = Argentina,
                   markers = spec_markers,
                   mutation = list(model = "equal", rate = 1e-3))
  )
  res <- evaluate_evidence(case)

  ref <- lr_distribution(case$genetic$models)
  expect_equal(res$combined_lr_dist$log10_lr, ref$log10_lr,
               tolerance = 1e-12)
  expect_equal(res$combined_lr_dist$p_h1, ref$p_h1, tolerance = 1e-12)
  expect_equal(res$combined_lr_dist$p_h2, ref$p_h2, tolerance = 1e-12)

  expect_equal(res$expected_log10_lr_h1,
               sum(res$per_feature_kl$e_log10_lr_h1),
               tolerance = 1e-9)
})

test_that("evaluate_evidence wires the thread policy", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  case <- mp_case(genetic = list(ped = ped, freqs = Argentina,
                                 markers = "D3S1358",
                                 mutation = list(model = "equal",
                                                 rate = 1e-3)))
  res <- evaluate_evidence(case, n_threads = 1L)
  expect_equal(res$n_threads, 1L)
})

test_that("non-genetic-only case evaluates and matches per-feature kernel", {
  feats <- list(sex = sex_feature("F", 0.05),
                sex2 = sex_feature("M", 0.10))
  case <- mp_case(non_genetic = feats)
  res <- evaluate_evidence(case)

  expect_equal(nrow(res$per_feature_kl), 2L)
  expect_true(all(res$per_feature_kl$kind == "non_genetic"))
  k1 <- per_feature_kl_cpp_wrap(feats[[1]])
  expect_equal(res$per_feature_kl$e_log10_lr_h1[1], k1$e_log10_lr_h1,
               tolerance = 1e-12)
})

test_that("print methods are stable", {
  skip_if_no_pedtools()
  ped <- pedtools::nuclearPed(1)
  case <- mp_case(
    genetic = list(ped = ped, freqs = Argentina, markers = "D3S1358",
                   mutation = list(model = "equal", rate = 1e-3)),
    non_genetic = list(sex = sex_feature())
  )
  expect_output(print(case), "<mp_case>")
  res <- evaluate_evidence(case)
  expect_output(print(res), "<evidence_result>")
  expect_output(summary(res), "Evidence evaluation summary")
})
