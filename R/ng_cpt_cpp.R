## F6.3 — R-side wrapper for the C++ kernel nongenetic_cpt_cpp().
##
## Mirrors ng_cpt_R() (R/r_ref_nongenetic.R): consumes a
## `nongenetic_feature`, flattens it into the POD inputs the C++ engine
## expects, runs the kernel, and rebuilds the labelled data.frame
## (`state`, `p_h1`, `p_h2`) with the same attributes so the F6.3 / F6.6
## cross-checks compare element-by-element against the R reference.

#' @noRd
ng_cpt_cpp_wrap <- function(feature) {
  if (!inherits(feature, "nongenetic_feature")) {
    stop("`feature` must be a 'nongenetic_feature' object.", call. = FALSE)
  }
  fclass <- feature$feature_class

  ## POD defaults (only the class-relevant fields are read by the core).
  fc <- 0L
  n_categories <- 0L
  error_is_matrix <- FALSE
  error_matrix <- numeric(0)
  error_scalar <- 0
  observed_index <- 0L
  reference_uniform <- FALSE
  reference_freqs <- numeric(0)
  range_lo <- 0
  range_hi <- 0
  sample_vec <- numeric(0)
  observed_value <- 0
  n_bins <- 0L
  alpha <- numeric(0)
  search_open <- TRUE

  if (fclass == "categorical") {
    fc <- 0L
    cats <- feature$categories
    n_categories <- length(cats)
    if (is.matrix(feature$error)) {
      error_is_matrix <- TRUE
      ## Positional K x K, row-major flatten (R matrices are
      ## column-major; transpose before vectorising).
      error_matrix <- as.numeric(t(unname(feature$error)))
    } else {
      error_scalar <- as.numeric(feature$error)
    }
    t_idx <- match(as.character(feature$observed), cats)
    if (is.na(t_idx)) {
      stop("`observed` is not one of the feature categories.",
           call. = FALSE)
    }
    observed_index <- as.integer(t_idx - 1L)
    reference_uniform <- identical(feature$model$reference, "uniform")
    reference_freqs <- as.numeric(feature$db_or_freqs)
    state <- cats
  } else if (fclass == "continuous") {
    fc <- 1L
    reference_uniform <- identical(feature$model$reference, "uniform")
    if (reference_uniform) {
      rg <- feature$model$range
      range_lo <- as.numeric(rg[1L])
      range_hi <- as.numeric(rg[2L])
    } else {
      sample_vec <- as.numeric(feature$db_or_freqs)
    }
    error_scalar <- as.numeric(feature$error)
    observed_value <- as.numeric(feature$observed)
    state <- NULL  # filled from the returned grid
  } else {
    fc <- 2L
    cuts <- feature$model$cuts
    n_bins <- length(cuts) + 1L
    alpha <- as.numeric(feature$error)
    search_open <- identical(feature$model$search, "open")
    if (!search_open) {
      db <- feature$db_or_freqs
      if (is.data.frame(db) || !is.numeric(db)) {
        stop("the non-genetic reference engine requires pre-binned ",
             "non-negative bin frequencies (length = length(cuts) + 1) ",
             "for a closed date search; the legacy stochastic ",
             "data.frame path stays in lr_birthdate().", call. = FALSE)
      }
      reference_freqs <- as.numeric(db)
    }
    state <- ng_date_bin_labels(cuts)
  }

  res <- nongenetic_cpt_cpp(
    feature_class = fc,
    n_categories = n_categories,
    error_is_matrix = error_is_matrix,
    error_matrix = error_matrix,
    error_scalar = error_scalar,
    observed_index = observed_index,
    reference_uniform = reference_uniform,
    reference_freqs = reference_freqs,
    range_lo = range_lo,
    range_hi = range_hi,
    sample = sample_vec,
    observed_value = observed_value,
    n_bins = n_bins,
    alpha = alpha,
    search_open = search_open
  )

  if (fclass == "continuous") {
    state <- as.character(res$grid)
  }

  out <- data.frame(state = state, p_h1 = res$p_h1, p_h2 = res$p_h2,
                     stringsAsFactors = FALSE)
  attr(out, "feature_type") <- feature$type
  attr(out, "observed") <- feature$observed
  out
}
