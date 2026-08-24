## F6.3 / F6.4 — R-side wrappers for the C++ non-genetic kernels.
##
## ng_cpt_cpp_wrap()            -> nongenetic_cpt_cpp()                 (F6.3)
## per_feature_kl_cpp_wrap()    -> nongenetic_per_feature_kl_cpp()      (F6.4)
## per_feature_lr_dist_cpp_wrap -> nongenetic_per_feature_lr_dist_cpp() (F6.4)
##
## All three flatten a `nongenetic_feature` into the identical POD input
## list via ng_feature_to_pod() and rebuild the labelled R object so the
## F6.3 / F6.4 / F6.6 cross-checks compare element-by-element against the
## R reference engine (ng_cpt_R / per_feature_kl_R / per_feature_lr_dist_R
## in r_ref_nongenetic.R).

## Flatten a `nongenetic_feature` into the POD argument list the C++
## kernels expect, plus the state labels for the categorical / date
## classes (continuous labels come back from the kernel grid). The
## class-irrelevant fields keep neutral defaults; only the fields read by
## the active feature class matter to the core.
#' @noRd
ng_feature_to_pod <- function(feature) {
  if (!inherits(feature, "nongenetic_feature")) {
    stop("`feature` must be a 'nongenetic_feature' object.", call. = FALSE)
  }
  fclass <- feature$feature_class

  args <- list(
    feature_class     = 0L,
    n_categories      = 0L,
    error_is_matrix   = FALSE,
    error_matrix      = numeric(0),
    error_scalar      = 0,
    observed_index    = 0L,
    reference_uniform = FALSE,
    reference_freqs   = numeric(0),
    range_lo          = 0,
    range_hi          = 0,
    sample            = numeric(0),
    observed_value    = 0,
    n_bins            = 0L,
    alpha             = numeric(0),
    search_open       = TRUE
  )
  state <- NULL

  if (fclass == "categorical") {
    args$feature_class <- 0L
    cats <- feature$categories
    args$n_categories <- length(cats)
    if (is.matrix(feature$error)) {
      args$error_is_matrix <- TRUE
      ## Positional K x K, row-major flatten (R matrices are
      ## column-major; transpose before vectorising).
      args$error_matrix <- as.numeric(t(unname(feature$error)))
    } else {
      args$error_scalar <- as.numeric(feature$error)
    }
    t_idx <- match(as.character(feature$observed), cats)
    if (is.na(t_idx)) {
      stop("`observed` is not one of the feature categories.",
           call. = FALSE)
    }
    args$observed_index <- as.integer(t_idx - 1L)
    args$reference_uniform <- identical(feature$model$reference, "uniform")
    args$reference_freqs <- as.numeric(feature$db_or_freqs)
    state <- cats
  } else if (fclass == "continuous") {
    args$feature_class <- 1L
    args$reference_uniform <- identical(feature$model$reference, "uniform")
    if (args$reference_uniform) {
      rg <- feature$model$range
      args$range_lo <- as.numeric(rg[1L])
      args$range_hi <- as.numeric(rg[2L])
    } else {
      args$sample <- as.numeric(feature$db_or_freqs)
    }
    args$error_scalar <- as.numeric(feature$error)
    args$observed_value <- as.numeric(feature$observed)
    state <- NULL  # filled from the returned grid
  } else {
    args$feature_class <- 2L
    cuts <- feature$model$cuts
    args$n_bins <- length(cuts) + 1L
    args$alpha <- as.numeric(feature$error)
    args$search_open <- identical(feature$model$search, "open")
    if (!args$search_open) {
      db <- feature$db_or_freqs
      if (is.data.frame(db) || !is.numeric(db)) {
        stop("the non-genetic reference engine requires pre-binned ",
             "non-negative bin frequencies (length = length(cuts) + 1) ",
             "for a closed date search; the legacy stochastic ",
             "data.frame path stays in lr_birthdate().", call. = FALSE)
      }
      args$reference_freqs <- as.numeric(db)
    }
    state <- ng_date_bin_labels(cuts)
  }

  list(args = args, state = state)
}

#' @noRd
ng_cpt_cpp_wrap <- function(feature) {
  pod <- ng_feature_to_pod(feature)
  res <- do.call(nongenetic_cpt_cpp, pod$args)

  state <- pod$state
  if (feature$feature_class == "continuous") {
    state <- as.character(res$grid)
  }

  out <- data.frame(state = state, p_h1 = res$p_h1, p_h2 = res$p_h2,
                     stringsAsFactors = FALSE)
  attr(out, "feature_type") <- feature$type
  attr(out, "observed") <- feature$observed
  out
}

## Bidirectional KL + expected log10 LR for one non-genetic feature.
## Mirrors per_feature_kl_R(): same four columns, `feature` = the
## feature type. The KLde-style absolute-continuity diagnostics returned
## by the kernel are not part of the R-ref schema and are dropped here
## (they remain available through the raw binding if ever needed).
#' @noRd
per_feature_kl_cpp_wrap <- function(feature) {
  pod <- ng_feature_to_pod(feature)
  r <- do.call(nongenetic_per_feature_kl_cpp, pod$args)
  data.frame(
    feature = feature$type,
    e_log10_lr_h1 = r$e_log10_lr_h1,
    e_log10_lr_h2 = r$e_log10_lr_h2,
    kl_h1_to_h2 = r$kl_h1_to_h2,
    kl_h2_to_h1 = r$kl_h2_to_h1,
    stringsAsFactors = FALSE
  )
}

## Sparse per-feature LR distribution. Mirrors per_feature_lr_dist_R():
## (log10_lr, p_h1, p_h2) with the feature_type / observed attributes.
#' @noRd
per_feature_lr_dist_cpp_wrap <- function(feature, aggregate = TRUE) {
  pod <- ng_feature_to_pod(feature)
  r <- do.call(nongenetic_per_feature_lr_dist_cpp,
               c(pod$args, list(aggregate = aggregate)))
  out <- data.frame(
    log10_lr = r$log10_lr,
    p_h1 = r$p_h1,
    p_h2 = r$p_h2,
    stringsAsFactors = FALSE
  )
  attr(out, "feature_type") <- feature$type
  attr(out, "observed") <- feature$observed
  out
}
