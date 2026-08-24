## R reference engine — per-feature non-genetic CPT, KL and LR
## distribution (F6.2). Internal-only helpers (no @export); the C++
## engine (F6.3--F6.5) is validated against these as oracles in F6.6.
## Slow / dense by design, like the genetic `r_ref_*` engines.
##
## This is the non-genetic generalisation requested by F6.2: the legacy
## functions `lr_sex`, `lr_age`, `lr_hair_color`, `lr_pigmentation`,
## `lr_birthdate`, `cpt_missing_person`, `cpt_population` all reduce to
## the same shape as the genetic per-marker path --- a discrete
## observation distribution under H1 (the missing person, observed with
## error) and under H2 (the population / reference marginal), with
## `LR(o) = P(o | H1) / P(o | H2)`. Downstream KL and LR distribution
## reuse the SAME generic helpers as the genetic engine
## (`log10_lr_from_probs`, `weighted_log10_lr_sum`, `aggregate_lr_dist`,
## defined in r_ref_per_marker.R), so the per-feature schema is
## bit-identical in shape to per-marker.
##
## Self-oracle reductions (exercised in F6.6):
##   categorical  sex / hair / eyes / region / pigmentation
##                K x K row-stochastic confusion E; `observed` is the
##                missing person's TRUE category (the `MPs` / `MPc`
##                argument of the legacy simulators). p_h1 = E[true, ],
##                p_h2 = population (or uniform). With K = 2 and a
##                scalar `eps` this is exactly `lr_sex`'s
##                (1 - eps) / Ps[MP] vs eps / Ps[other]; with a full
##                matrix it is exactly `lr_hair_color`'s
##                epc[MPc, o] / Pc[o].
##   continuous   age: a population grid (integer ages over `range`, or
##                the empirical sample) with a scalar mis-binning rate.
##                Collapsing the grid to the two cells {within tol,
##                outside} reproduces `lr_age`'s T1 / T0 binary LR
##                exactly (the verifier builds that 2-category feature
##                in F6.6).
##   date         birthdate: the deterministic Dirichlet *mean*
##                p_h1 = alpha / sum(alpha) over the signed
##                declared-minus-actual discrepancy bins. `lr_birthdate`
##                estimates the same quantity stochastically (method of
##                moments over `DirichletReg::rdirichlet` draws), which
##                converges to alpha / sum(alpha); the C++ engine in
##                F6.3 must be deterministic, so the reference is the
##                limit, not the sampler.
##
## Boundary convention is shared with the genetic engine: states with
## p_h1 > 0, p_h2 = 0 give log10 LR = +Inf; p_h1 = 0, p_h2 > 0 give
## -Inf; sums use the limit 0 * log(0 / x) = 0.

#' @noRd
ng_confusion_matrix <- function(error, K) {
  if (is.matrix(error)) {
    ## Positional: row i / column j are the i-th / j-th category in the
    ## order of `db_or_freqs` (validated in nongenetic_feature()).
    return(unname(error))
  }
  eps <- error
  if (K < 2L) {
    stop("a categorical feature needs at least two categories.",
         call. = FALSE)
  }
  m <- matrix(eps / (K - 1L), nrow = K, ncol = K)
  diag(m) <- 1 - eps
  m
}

#' @noRd
ng_cpt_categorical <- function(feature) {
  cats <- feature$categories
  K <- length(cats)
  E <- ng_confusion_matrix(feature$error, K)
  ## `observed` is stored as a character category label by the
  ## constructor; it plays the role of the missing person's TRUE
  ## category (the conditioning value for H1).
  t_idx <- match(as.character(feature$observed), cats)
  if (is.na(t_idx)) {
    stop("`observed` is not one of the feature categories.", call. = FALSE)
  }
  p_h1 <- as.numeric(E[t_idx, ])
  ref <- feature$model$reference
  p_h2 <- if (identical(ref, "uniform")) {
    rep(1 / K, K)
  } else {
    as.numeric(feature$db_or_freqs)
  }
  data.frame(state = cats, p_h1 = p_h1, p_h2 = p_h2,
             stringsAsFactors = FALSE)
}

#' @noRd
ng_cpt_continuous <- function(feature) {
  ref <- feature$model$reference
  if (identical(ref, "uniform")) {
    rg <- feature$model$range
    grid <- seq(floor(rg[1L]), ceiling(rg[2L]))
    p_h2 <- rep(1 / length(grid), length(grid))
  } else {
    s <- feature$db_or_freqs
    grid <- sort(unique(s))
    counts <- tabulate(match(s, grid), nbins = length(grid))
    p_h2 <- counts / sum(counts)
  }
  G <- length(grid)
  eps <- feature$error
  ti <- which.min(abs(grid - feature$observed))
  if (G < 2L) {
    p_h1 <- 1
  } else {
    p_h1 <- rep(eps / (G - 1L), G)
    p_h1[ti] <- 1 - eps
  }
  data.frame(state = as.character(grid), p_h1 = p_h1, p_h2 = p_h2,
             stringsAsFactors = FALSE)
}

#' @noRd
ng_date_bin_labels <- function(cuts) {
  brks <- c(-Inf, cuts, Inf)
  paste0("(", utils::head(brks, -1L), ",", brks[-1L], "]")
}

#' @noRd
ng_cpt_date <- function(feature) {
  cuts <- feature$model$cuts
  nbins <- length(cuts) + 1L
  alpha <- feature$error
  p_h1 <- alpha / sum(alpha)
  if (identical(feature$model$search, "open")) {
    p_h2 <- rep(1 / nbins, nbins)
  } else {
    db <- feature$db_or_freqs
    if (is.data.frame(db) || !is.numeric(db)) {
      stop("the non-genetic reference engine requires pre-binned ",
           "non-negative bin frequencies (length = length(cuts) + 1) ",
           "for a closed date search; the legacy stochastic ",
           "data.frame path stays in lr_birthdate().", call. = FALSE)
    }
    if (sum(db) <= 0) {
      stop("closed-search bin frequencies must have positive total ",
           "mass.", call. = FALSE)
    }
    p_h2 <- db / sum(db)
  }
  data.frame(state = ng_date_bin_labels(cuts),
             p_h1 = as.numeric(p_h1), p_h2 = as.numeric(p_h2),
             stringsAsFactors = FALSE)
}

#' @noRd
ng_cpt_R <- function(feature) {
  if (!inherits(feature, "nongenetic_feature")) {
    stop("`feature` must be a 'nongenetic_feature' object.", call. = FALSE)
  }
  out <- switch(feature$feature_class,
    categorical = ng_cpt_categorical(feature),
    continuous  = ng_cpt_continuous(feature),
    date        = ng_cpt_date(feature)
  )
  attr(out, "feature_type") <- feature$type
  attr(out, "observed") <- feature$observed
  out
}

#' @noRd
per_feature_lr_dist_R <- function(feature, aggregate = TRUE) {
  cpt <- ng_cpt_R(feature)
  log10_lr <- log10_lr_from_probs(cpt$p_h1, cpt$p_h2)
  out <- data.frame(
    log10_lr = log10_lr,
    p_h1 = cpt$p_h1,
    p_h2 = cpt$p_h2
  )
  if (aggregate) {
    out <- aggregate_lr_dist(out)
  }
  attr(out, "feature_type") <- attr(cpt, "feature_type")
  attr(out, "observed") <- attr(cpt, "observed")
  out
}

#' @noRd
per_feature_kl_R <- function(feature) {
  cpt <- ng_cpt_R(feature)
  P1 <- cpt$p_h1
  P2 <- cpt$p_h2
  log10_lr <- log10_lr_from_probs(P1, P2)
  ln10 <- log(10)
  e_h1 <- weighted_log10_lr_sum(P1, log10_lr)
  e_h2 <- weighted_log10_lr_sum(P2, log10_lr)
  data.frame(
    feature = attr(cpt, "feature_type"),
    e_log10_lr_h1 = e_h1,
    e_log10_lr_h2 = e_h2,
    kl_h1_to_h2 = e_h1 * ln10,
    kl_h2_to_h1 = -e_h2 * ln10,
    stringsAsFactors = FALSE
  )
}

## ---------------------------------------------------------------------
## Soft deprecation of the legacy non-genetic API (F6.2).
##
## The seven legacy entry points keep their exact behaviour (so the
## existing test suite, the vignettes, and the F6.6 self-oracle all stay
## green) but now generalise to `nongenetic_feature()` + the engine
## above. They emit a single, suppressible message per session pointing
## users at the unified constructor --- "soft" deprecation for the 2.0
## release-candidate cycle, not a hard `.Deprecated()` warning.
## ---------------------------------------------------------------------

.mispitools_dep_state <- new.env(parent = emptyenv())

#' @noRd
ng_reset_deprecation_state <- function() {
  rm(list = ls(envir = .mispitools_dep_state, all.names = TRUE),
     envir = .mispitools_dep_state)
  invisible()
}

#' @noRd
ng_soft_deprecate <- function(old, hint) {
  if (isTRUE(getOption("mispitools.suppress_deprecation", FALSE))) {
    return(invisible())
  }
  key <- paste0("warned_", old)
  if (isTRUE(.mispitools_dep_state[[key]])) {
    return(invisible())
  }
  .mispitools_dep_state[[key]] <- TRUE
  message(sprintf(
    paste0("%s() is soft-deprecated in mispitools 2.0 and will be ",
           "removed after the 2.0 release cycle. %s\n",
           "Suppress this notice with ",
           "options(mispitools.suppress_deprecation = TRUE)."),
    old, hint))
  invisible()
}
