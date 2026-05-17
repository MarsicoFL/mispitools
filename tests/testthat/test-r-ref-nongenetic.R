## F6.2 -- R-pure non-genetic per-feature reference engine and the
## soft-deprecation of the legacy non-genetic API. These engines are the
## oracle the C++ kernel (F6.3--F6.5) must replicate and the self-oracle
## the legacy functions are cross-checked against in F6.6.

ln10 <- log(10)

# --- ng_cpt_R: structural contract -----------------------------------

test_that("ng_cpt_R rejects non-features and returns the joint schema", {
  expect_error(ng_cpt_R(list()), "nongenetic_feature")

  f <- nongenetic_feature("sex", observed = "F",
                          db_or_freqs = c(F = 0.5, M = 0.5), error = 0.05)
  cpt <- ng_cpt_R(f)
  expect_s3_class(cpt, "data.frame")
  expect_named(cpt, c("state", "p_h1", "p_h2"))
  expect_equal(attr(cpt, "feature_type"), "sex")
  ## p_h1 is a proper conditional distribution; p_h2 a proper marginal.
  expect_equal(sum(cpt$p_h1), 1, tolerance = 1e-12)
  expect_equal(sum(cpt$p_h2), 1, tolerance = 1e-12)
})

# --- categorical: reproduces lr_sex exactly --------------------------

test_that("categorical engine reproduces lr_sex's deterministic LR", {
  eps <- 0.05
  Ps <- c(0.5, 0.5)
  f <- nongenetic_feature("sex", observed = "F",
                          db_or_freqs = c(F = Ps[1], M = Ps[2]),
                          error = eps)
  cpt <- ng_cpt_R(f)
  lr <- cpt$p_h1 / cpt$p_h2
  names(lr) <- cpt$state
  ## lr_sex: match -> (1 - eps) / Ps[MP]; other -> eps / Ps[other].
  expect_equal(unname(lr["F"]), (1 - eps) / Ps[1], tolerance = 1e-12)
  expect_equal(unname(lr["M"]), eps / Ps[2], tolerance = 1e-12)

  ## Uniform reference => H2 is 1/K, not the population.
  fu <- nongenetic_feature("sex", observed = "M",
                           model = list(reference = "uniform"),
                           db_or_freqs = c(F = 0.7, M = 0.3),
                           error = eps)
  cu <- ng_cpt_R(fu)
  expect_equal(cu$p_h2, c(0.5, 0.5), tolerance = 1e-12)
})

# --- categorical: reproduces lr_hair_color exactly -------------------

test_that("categorical engine reproduces lr_hair_color's epc[MPc,o]/Pc[o]", {
  E <- error_matrix_hair()
  Pc <- c(0.3, 0.2, 0.25, 0.15, 0.1)
  for (mpc in 1:5) {
    f <- nongenetic_feature(
      "hair", observed = mpc,
      db_or_freqs = stats::setNames(Pc, as.character(1:5)),
      error = E)
    cpt <- ng_cpt_R(f)
    expect_equal(cpt$p_h1, unname(E[mpc, ]), tolerance = 1e-12)
    expect_equal(cpt$p_h1 / cpt$p_h2, unname(E[mpc, ]) / Pc,
                 tolerance = 1e-12)
  }
})

# --- continuous (age) -------------------------------------------------

test_that("continuous engine puts (1-eps) on the observed bin, uniform H2", {
  eps <- 0.05
  f <- nongenetic_feature("age", observed = 42,
                          model = list(reference = "uniform",
                                       range = c(1, 80)),
                          error = eps)
  cpt <- ng_cpt_R(f)
  expect_equal(nrow(cpt), 80L)
  expect_equal(sum(cpt$p_h1), 1, tolerance = 1e-12)
  expect_equal(sum(cpt$p_h2), 1, tolerance = 1e-12)
  ti <- which(cpt$state == "42")
  expect_equal(cpt$p_h1[ti], 1 - eps, tolerance = 1e-12)
  expect_true(all(abs(cpt$p_h1[-ti] - eps / 79) < 1e-12))
  expect_equal(cpt$p_h2, rep(1 / 80, 80), tolerance = 1e-12)
})

test_that("continuous engine 2-cell collapse equals lr_age's T1/T0 LR", {
  ## lr_age(MPa = 40, MPr = 6, epa = 0.05, modelA = "uniform"):
  ##   P(T1) = 2*MPr/80, LR1 = (1-epa)/P(T1), LR0 = epa/P(T0).
  ## The unified analogue is a 2-category feature {T1, T0}.
  epa <- 0.05
  T1p <- (46 - 34) / 80
  f <- nongenetic_feature(
    type = "custom", observed = "T1",
    model = list(class = "categorical", reference = "marginal"),
    db_or_freqs = c(T1 = T1p, T0 = 1 - T1p),
    error = epa)
  cpt <- ng_cpt_R(f)
  lr <- cpt$p_h1 / cpt$p_h2
  names(lr) <- cpt$state
  expect_equal(unname(lr["T1"]), (1 - epa) / T1p, tolerance = 1e-12)
  expect_equal(unname(lr["T0"]), epa / (1 - T1p), tolerance = 1e-12)
})

# --- date (birthdate) -------------------------------------------------

test_that("date engine is the deterministic Dirichlet mean, open => uniform H2", {
  alpha <- c(1, 4, 60, 11, 6, 4, 4)
  f <- nongenetic_feature("birthdate", observed = 45, error = alpha)
  cpt <- ng_cpt_R(f)
  expect_equal(nrow(cpt), length(alpha))
  expect_equal(cpt$p_h1, alpha / sum(alpha), tolerance = 1e-12)
  expect_equal(cpt$p_h2, rep(1 / length(alpha), length(alpha)),
               tolerance = 1e-12)
  expect_equal(cpt$p_h1 / cpt$p_h2,
               (alpha / sum(alpha)) * length(alpha), tolerance = 1e-12)
})

test_that("date engine closed search uses normalised bin frequencies", {
  alpha <- c(1, 4, 60, 11, 6, 4, 4)
  db <- c(2, 3, 50, 10, 8, 5, 2)
  f <- nongenetic_feature(
    "birthdate", observed = 45,
    model = list(search = "closed",
                 cuts = c(-120, -30, 30, 120, 240, 360)),
    db_or_freqs = db, error = alpha)
  cpt <- ng_cpt_R(f)
  expect_equal(cpt$p_h2, db / sum(db), tolerance = 1e-12)
})

# --- downstream KL / LR distribution symmetry ------------------------

test_that("per_feature_kl_R is self-consistent and matches per-marker schema", {
  f <- nongenetic_feature("hair", observed = 2,
                          db_or_freqs = c("1" = 0.3, "2" = 0.2, "3" = 0.25,
                                          "4" = 0.15, "5" = 0.1),
                          error = error_matrix_hair())
  kl <- per_feature_kl_R(f)
  expect_named(kl, c("feature", "e_log10_lr_h1", "e_log10_lr_h2",
                     "kl_h1_to_h2", "kl_h2_to_h1"))
  expect_equal(kl$kl_h1_to_h2, kl$e_log10_lr_h1 * ln10, tolerance = 1e-12)
  expect_equal(kl$kl_h2_to_h1, -kl$e_log10_lr_h2 * ln10, tolerance = 1e-12)
  ## KL(P1||P2) and KL(P2||P1) are non-negative for these strictly
  ## positive distributions.
  expect_gte(kl$kl_h1_to_h2, -1e-12)
  expect_gte(kl$kl_h2_to_h1, -1e-12)
})

test_that("per_feature_lr_dist_R is sorted, aggregated and sums to 1", {
  f <- nongenetic_feature("sex", observed = "F",
                          db_or_freqs = c(F = 0.5, M = 0.5), error = 0.05)
  d <- per_feature_lr_dist_R(f)
  expect_named(d, c("log10_lr", "p_h1", "p_h2"))
  expect_false(is.unsorted(d$log10_lr))
  expect_equal(sum(d$p_h1), 1, tolerance = 1e-12)
  expect_equal(sum(d$p_h2), 1, tolerance = 1e-12)
  ## Mean log10 LR under H1 must equal the KL-engine expectation.
  m <- sum(d$p_h1 * d$log10_lr)
  expect_equal(m, per_feature_kl_R(f)$e_log10_lr_h1, tolerance = 1e-12)
})

# --- soft deprecation of the legacy API ------------------------------

test_that("ng_soft_deprecate fires once per session and is suppressible", {
  op0 <- options(mispitools.suppress_deprecation = FALSE)
  on.exit(options(op0), add = TRUE)
  ng_reset_deprecation_state()
  expect_message(ng_soft_deprecate("demo_fn", "Use new()."),
                 "soft-deprecated")
  ## Second call in the same session is silent.
  expect_message(ng_soft_deprecate("demo_fn", "Use new()."), NA)

  ng_reset_deprecation_state()
  op <- options(mispitools.suppress_deprecation = TRUE)
  on.exit(options(op), add = TRUE)
  expect_message(ng_soft_deprecate("demo_fn", "Use new()."), NA)
})

test_that("legacy lr_sex stays numerically intact and only soft-warns", {
  op0 <- options(mispitools.suppress_deprecation = FALSE)
  on.exit(options(op0), add = TRUE)
  ng_reset_deprecation_state()
  expect_message(
    res <- lr_sex(MPs = "F", eps = 0.05, Ps = c(0.5, 0.5),
                  numsims = 50, LR = TRUE, seed = 1234),
    "soft-deprecated")
  expect_s3_class(res, "data.frame")
  expect_true(all(c("Sexo", "LRs") %in% names(res)))
  ## Behaviour unchanged: matching obs -> (1-eps)/Ps, else eps/Ps.
  expect_true(all(res$LRs[res$Sexo == "F"] == (1 - 0.05) / 0.5))
  expect_true(all(res$LRs[res$Sexo == "M"] == 0.05 / 0.5))
  ## The engine reproduces those same two LR values.
  cpt <- ng_cpt_R(nongenetic_feature("sex", observed = "F",
                                     db_or_freqs = c(F = 0.5, M = 0.5),
                                     error = 0.05))
  expect_equal(sort(unique(res$LRs)),
               sort((cpt$p_h1 / cpt$p_h2)), tolerance = 1e-12)
})
