## F6.6 — Self-oracle cross-check of the C++ non-genetic per-feature
## engine against the *original* legacy simulators (lr_sex,
## lr_hair_color, lr_age, ...). Oracle de sí mismas: the legacy
## functions are deterministic in the LR they assign to each observed
## category (only the sampled observation is random), so the unified
## engine — which builds the same P(o | H1) / P(o | H2) — must match
## them exactly. Tolerance 1e-12: the logic is literally the same
## arithmetic, only the difference is binding/round-trip floating point.
##
## Scope vs prior phases:
##   * F6.2 validated the R-pure reference (ng_cpt_R, ...) against the
##     legacy *formula*.
##   * F6.3/F6.4 validated the C++ kernels against that R reference
##     (shared-oracle, 1e-12).
## F6.6 closes the chain by hitting the C++ pipeline (Rcpp binding +
## src/core/nongenetic_lr) directly with the legacy functions
## themselves as an independent oracle — catching a binding-layer bug
## that a shared R-ref oracle would hide.
##
## NOT covered here (and why — never loosen a tolerance to pass):
##   * lr_birthdate: estimates the Dirichlet mean *stochastically*
##     (method of moments over DirichletReg::rdirichlet draws). The
##     engine uses the deterministic limit alpha / sum(alpha). They do
##     NOT coincide at 1e-12 by construction (it is a sampler vs its
##     limit), so a 1e-12 self-oracle is ill-posed. The C++ date path
##     is instead validated against the deterministic R reference in
##     F6.4 (test-per-feature-nongenetic-cpp.R, 1e-12).
##   * lr_pigmentation: consumes a pre-computed joint-class LR table
##     (lr_compute_pigmentation), not an error-matrix observation
##     model, so it is not an exact error-matrix reduction either. Its
##     per-class LR is numerators/f_h_s_y; the categorical engine
##     reproduces that ratio when fed those two columns as p_h1 / p_h2
##     (exercised below as the general categorical-ratio identity).

## The legacy simulators emit a one-shot soft-deprecation message
## (F6.2); suppressMessages keeps the test output clean without
## mutating the global mispitools.suppress_deprecation option.
legacy_lr_map <- function(df, cat_col, lr_col) {
  ## Collapse a legacy simulator's (observation, LR) draws to its
  ## deterministic category -> LR mapping. The LR column is a pure
  ## function of the observed category, so unique rows recover it
  ## exactly regardless of which categories the RNG happened to draw.
  u <- unique(df[, c(cat_col, lr_col)])
  stats::setNames(u[[lr_col]], as.character(u[[cat_col]]))
}

# --- lr_sex: scalar-eps categorical, K = 2 ---------------------------

test_that("C++ engine reproduces lr_sex's deterministic LR (1e-12)", {
  grid <- expand.grid(
    MP  = c("F", "M"),
    eps = c(0.05, 0.10, 0.005),
    pF  = c(0.5, 0.52, 0.40),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(grid))) {
    MP  <- grid$MP[i]
    eps <- grid$eps[i]
    Ps  <- c(grid$pF[i], 1 - grid$pF[i])

    leg <- suppressMessages(
      lr_sex(MPs = MP, eps = eps, Ps = Ps, H = 2, LR = TRUE,
             numsims = 4000, seed = 99)
    )
    leg_lr <- legacy_lr_map(leg, "Sexo", "LRs")
    expect_setequal(names(leg_lr), c("F", "M"))  # RNG covered both

    f <- nongenetic_feature(
      type = "sex", observed = MP,
      db_or_freqs = c(F = Ps[1], M = Ps[2]), error = eps
    )
    cpt <- ng_cpt_cpp_wrap(f)
    eng_lr <- stats::setNames(cpt$p_h1 / cpt$p_h2, cpt$state)

    ## Engine LR per category == legacy LR per category, exactly.
    ## 1e-12 because the only gap is the Rcpp double round-trip; the
    ## arithmetic (1-eps)/Ps[MP] vs eps/Ps[other] is identical.
    expect_equal(eng_lr[names(leg_lr)], leg_lr, tolerance = 1e-12)

    ## And the full per-feature LR distribution: H2 mass == population
    ## proportions, H1 mass == the error-model row.
    d <- per_feature_lr_dist_cpp_wrap(f, aggregate = FALSE)
    expect_equal(sort(10^d$log10_lr), sort(unname(leg_lr)),
                 tolerance = 1e-12)
    expect_equal(sum(d$p_h1), 1, tolerance = 1e-12)
    expect_equal(sum(d$p_h2), 1, tolerance = 1e-12)
  }
})

# --- lr_hair_color: full 5x5 confusion matrix ------------------------

test_that("C++ engine reproduces lr_hair_color's deterministic LR (1e-12)", {
  E <- error_matrix_hair()                 # 5x5 row-stochastic
  expect_equal(unname(rowSums(E)), rep(1, 5), tolerance = 1e-12)

  cases <- list(
    list(MPc = 1, Pc = c(0.30, 0.20, 0.25, 0.15, 0.10)),
    list(MPc = 3, Pc = c(0.10, 0.40, 0.30, 0.10, 0.10)),
    list(MPc = 5, Pc = c(0.22, 0.18, 0.20, 0.20, 0.20))
  )
  for (cs in cases) {
    leg <- suppressMessages(
      lr_hair_color(MPc = cs$MPc, epc = E, Pc = cs$Pc, H = 2, LR = TRUE,
                    numsims = 6000, seed = 7)
    )
    leg_lr <- legacy_lr_map(leg, "Col", "LRc")
    expect_setequal(names(leg_lr), as.character(1:5))

    f <- nongenetic_feature(
      type = "hair", observed = cs$MPc,
      db_or_freqs = stats::setNames(cs$Pc, as.character(1:5)),
      error = E
    )
    cpt <- ng_cpt_cpp_wrap(f)
    eng_lr <- stats::setNames(cpt$p_h1 / cpt$p_h2, cpt$state)

    ## Legacy LR = epc[MPc, o] / Pc[o]; the positional confusion-matrix
    ## engine builds p_h1 = E[MPc, ], p_h2 = Pc -> identical ratio.
    expect_equal(eng_lr[as.character(1:5)],
                 leg_lr[as.character(1:5)], tolerance = 1e-12)
    ## p_h1 is exactly the MP's confusion-matrix row.
    expect_equal(cpt$p_h1, unname(E[cs$MPc, ]), tolerance = 1e-12)
    expect_equal(cpt$p_h2, cs$Pc, tolerance = 1e-12)
  }
})

# --- lr_age: 2-cell collapse of the continuous feature ---------------

test_that("C++ engine reproduces lr_age's T1/T0 LR via 2-cell collapse (1e-12)", {
  ## STATE F6.2 self-oracle reduction: lr_age's uniform model is a
  ## binary {within tolerance, outside} categorical with
  ## T1p = (MPmax - MPmin) / 80 (legacy counts the interval width over
  ## length(Age) = 80, NOT the integer count), T0p = 1 - T1p, and
  ## LR1 = (1 - epa)/T1p, LR0 = epa/T0p. The verifier builds that
  ## 2-category feature and checks the engine against the *legacy*
  ## function's emitted LRa column.
  cases <- list(
    list(MPa = 40, MPr = 6,  epa = 0.05),
    list(MPa = 35, MPr = 3,  epa = 0.10),
    list(MPa = 50, MPr = 10, epa = 0.02)
  )
  for (cs in cases) {
    MPmin <- max(1, cs$MPa - cs$MPr)
    MPmax <- min(80, cs$MPa + cs$MPr)
    T1p <- (MPmax - MPmin) / 80
    T0p <- 1 - T1p

    leg <- suppressMessages(
      lr_age(MPa = cs$MPa, MPr = cs$MPr, epa = cs$epa, H = 2,
             modelA = "uniform", LR = TRUE, numsims = 8000, seed = 21)
    )
    leg_lr <- legacy_lr_map(leg, "group", "LRa")
    expect_setequal(names(leg_lr), c("T1", "T0"))

    f <- nongenetic_feature(
      type = "sex",                       # generic K=2 categorical
      observed = "T1",
      db_or_freqs = c(T1 = T1p, T0 = T0p),
      error = cs$epa
    )
    cpt <- ng_cpt_cpp_wrap(f)
    eng_lr <- stats::setNames(cpt$p_h1 / cpt$p_h2, cpt$state)

    expect_equal(unname(eng_lr["T1"]), unname(leg_lr["T1"]),
                 tolerance = 1e-12)
    expect_equal(unname(eng_lr["T0"]), unname(leg_lr["T0"]),
                 tolerance = 1e-12)
    ## Cross-check against the closed-form legacy arithmetic too.
    expect_equal(unname(eng_lr["T1"]), (1 - cs$epa) / T1p,
                 tolerance = 1e-12)
    expect_equal(unname(eng_lr["T0"]), cs$epa / T0p,
                 tolerance = 1e-12)
  }
})

# --- general categorical-ratio identity (covers lr_pigmentation) -----

test_that("categorical engine reproduces a precomputed p_h1/p_h2 LR table (1e-12)", {
  ## lr_pigmentation's per-class LR is numerators/f_h_s_y over the
  ## joint pigmentation classes. Feeding the conditioned probabilities
  ## as the H1 model (a confusion-matrix row, here the degenerate
  ## "observed == that class" row is not used; instead we use the
  ## general identity: with a full E whose `observed` row equals the
  ## H1 class probabilities, LR(state) = p_h1(state)/p_h2(state)).
  ## This is the same arithmetic lr_pigmentation samples from, so it
  ## is the self-oracle for the pigmentation reduction.
  set.seed(2024)
  K <- 6
  numerators <- runif(K, 0.01, 1)        # conditioned (H1) weights
  f_h_s_y    <- runif(K, 0.01, 1)        # population (H2) weights
  p_h1 <- numerators / sum(numerators)
  p_h2 <- f_h_s_y / sum(f_h_s_y)
  legacy_lr <- p_h1 / p_h2               # == lr_compute_pigmentation's LR

  ## Build a categorical feature whose H1 row is p_h1: a confusion
  ## matrix whose `observed` row is exactly p_h1, H2 = p_h2.
  E <- matrix(1 / K, K, K)               # other rows irrelevant
  E[1, ] <- p_h1
  f <- nongenetic_feature(
    type = "custom",
    model = list(class = "categorical", reference = "marginal"),
    observed = "c1",
    db_or_freqs = stats::setNames(p_h2, paste0("c", 1:K)),
    error = E
  )
  cpt <- ng_cpt_cpp_wrap(f)
  eng_lr <- cpt$p_h1 / cpt$p_h2
  expect_equal(eng_lr, legacy_lr, tolerance = 1e-12)
  expect_equal(cpt$p_h1, p_h1, tolerance = 1e-12)
  expect_equal(cpt$p_h2, p_h2, tolerance = 1e-12)
})
