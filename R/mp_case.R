
#' Build a unified missing-person evidence case
#'
#' @description
#' Bundles the genetic and non-genetic evidence for a single
#' missing-person comparison into one object that
#' [evaluate_evidence()] consumes. The genetic side is described by a
#' pedigree, a per-marker frequency source and a mutation / linkage
#' specification shared across markers; the non-genetic side is a list
#' of [nongenetic_feature()] objects. `combine` selects how the
#' per-feature likelihood ratios are pooled into the combined
#' distribution.
#'
#' @param genetic Either `NULL` (no genetic evidence) or a named list
#'   describing the genetic side:
#'   \describe{
#'     \item{`ped`}{A `pedtools::ped` object (required).}
#'     \item{`freqs`}{Allele frequencies (required). Either a named
#'       list of named numeric vectors (one per marker, the
#'       `pedtools` format produced by [get_allele_freqs()]) or a
#'       wide `data.frame` with a first `Allele` column followed by
#'       one column per marker (the built-in databases such as
#'       [Argentina]). Alleles with zero or missing frequency are
#'       dropped and each marker is renormalised to sum to 1.}
#'     \item{`mutation`}{Optional mutation model applied to every
#'       marker; see [marker_model()]. Defaults to
#'       `list(model = "none", rate = 0)`.}
#'     \item{`linkage`}{Optional linkage descriptor applied to every
#'       marker; see [marker_model()]. Defaults to `NULL`.}
#'     \item{`markers`}{Optional character vector selecting a subset of
#'       markers (by name) from `freqs`. `NULL` (default) uses every
#'       available marker.}
#'   }
#' @param non_genetic Either `NULL` (no non-genetic evidence) or a
#'   list of [nongenetic_feature()] objects. A single
#'   `nongenetic_feature` is accepted and wrapped in a list. List
#'   names, when present, label the features in the output.
#' @param combine Evidence-combination mode. `"independent"` (default)
#'   sums the per-feature `log10` LRs assuming conditional
#'   independence; `"markov_se"` follows the Egeland-Marsico (2026)
#'   first-order Markov chain over feature states and requires
#'   `transition`.
#' @param transition Required only when `combine = "markov_se"`: a
#'   list of `length(features) - 1` row-stochastic transition matrices
#'   in feature order. Ignored for `"independent"`.
#'
#' @return An object of class `"mp_case"`: a list with components
#'   `genetic` (with the built `models` list), `non_genetic` (the
#'   feature list), `combine`, `transition`, `n_markers` and
#'   `n_features`.
#'
#' @details
#' At least one of `genetic` or `non_genetic` must supply evidence.
#' Each marker frequency vector is cleaned (zero / `NA` alleles
#' dropped) and renormalised before a [marker_model()] is built, so
#' the wide population databases can be passed verbatim. Construction
#' is purely validating; no CPT is computed until
#' [evaluate_evidence()].
#'
#' @seealso [evaluate_evidence()], [marker_model()],
#'   [nongenetic_feature()], [get_allele_freqs()].
#'
#' @export
#' @examples
#' if (requireNamespace("pedtools", quietly = TRUE)) {
#'   ped <- pedtools::nuclearPed(1)
#'   case <- mp_case(
#'     genetic = list(
#'       ped = ped,
#'       freqs = Argentina,
#'       markers = c("D3S1358", "THO1"),
#'       mutation = list(model = "equal", rate = 1e-3)
#'     ),
#'     non_genetic = list(
#'       sex = nongenetic_feature(
#'         type = "sex", observed = "F",
#'         db_or_freqs = c(F = 0.5, M = 0.5), error = 0.05
#'       )
#'     )
#'   )
#'   print(case)
#' }
mp_case <- function(genetic = NULL,
                    non_genetic = NULL,
                    combine = c("independent", "markov_se"),
                    transition = NULL) {
  combine <- match.arg(combine)

  models <- list()
  genetic_spec <- NULL
  if (!is.null(genetic)) {
    genetic_spec <- validate_genetic_spec(genetic)
    models <- build_marker_models(genetic_spec)
  }

  features <- validate_non_genetic(non_genetic)

  if (length(models) == 0L && length(features) == 0L) {
    stop("`mp_case()` needs at least one source of evidence: supply ",
         "`genetic` and/or `non_genetic`.", call. = FALSE)
  }

  if (combine == "markov_se") {
    n_feat <- length(models) + length(features)
    if (is.null(transition) || !is.list(transition) ||
        length(transition) != max(n_feat - 1L, 0L)) {
      stop("`combine = \"markov_se\"` requires `transition`: a list of ",
           "length(features) - 1 row-stochastic matrices (got ",
           if (is.null(transition)) "NULL" else
             paste0("length ", length(transition)),
           ", expected ", max(n_feat - 1L, 0L), ").", call. = FALSE)
    }
  }

  structure(
    list(
      genetic = if (is.null(genetic_spec)) NULL else
        c(genetic_spec, list(models = models)),
      non_genetic = features,
      combine = combine,
      transition = transition,
      n_markers = length(models),
      n_features = length(features)
    ),
    class = "mp_case"
  )
}

#' @noRd
validate_genetic_spec <- function(genetic) {
  if (!is.list(genetic) || is.data.frame(genetic)) {
    stop("`genetic` must be a named list (ped, freqs, ...) or NULL.",
         call. = FALSE)
  }
  if (is.null(genetic$ped)) {
    stop("`genetic$ped` is required.", call. = FALSE)
  }
  if (is.null(genetic$freqs)) {
    stop("`genetic$freqs` is required.", call. = FALSE)
  }
  mutation <- genetic$mutation %||% list(model = "none", rate = 0)
  linkage <- genetic$linkage
  markers <- genetic$markers

  freqs_list <- if (is.data.frame(genetic$freqs)) {
    get_allele_freqs(genetic$freqs)
  } else if (is.list(genetic$freqs)) {
    genetic$freqs
  } else {
    stop("`genetic$freqs` must be a wide data.frame (Allele + marker ",
         "columns) or a named list of allele-frequency vectors.",
         call. = FALSE)
  }
  if (is.null(names(freqs_list)) || any(!nzchar(names(freqs_list)))) {
    stop("`genetic$freqs` must name every marker.", call. = FALSE)
  }

  if (!is.null(markers)) {
    if (!is.character(markers) || anyNA(markers) || !length(markers)) {
      stop("`genetic$markers` must be a non-empty character vector or ",
           "NULL.", call. = FALSE)
    }
    miss <- setdiff(markers, names(freqs_list))
    if (length(miss)) {
      stop("markers absent from `genetic$freqs`: ",
           paste(miss, collapse = ", "), ".", call. = FALSE)
    }
    freqs_list <- freqs_list[markers]
  }

  list(ped = genetic$ped, freqs = freqs_list, mutation = mutation,
       linkage = linkage, markers = names(freqs_list))
}

#' @noRd
clean_marker_freqs <- function(fr, marker) {
  fr <- unlist(fr, use.names = TRUE)
  if (is.null(names(fr)) || any(!nzchar(names(fr)))) {
    stop("frequencies for marker '", marker,
         "' must be named by allele.", call. = FALSE)
  }
  fr <- suppressWarnings(stats::setNames(as.numeric(fr), names(fr)))
  keep <- !is.na(fr) & fr > 0
  fr <- fr[keep]
  if (length(fr) < 2L) {
    stop("marker '", marker, "' has fewer than two alleles with ",
         "positive frequency.", call. = FALSE)
  }
  fr / sum(fr)
}

#' @noRd
build_marker_models <- function(spec) {
  ids <- spec$markers
  models <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    fr <- clean_marker_freqs(spec$freqs[[ids[i]]], ids[i])
    models[[i]] <- marker_model(
      ped = spec$ped,
      marker_id = ids[i],
      freqs = fr,
      mutation = spec$mutation,
      linkage = spec$linkage
    )
  }
  names(models) <- ids
  models
}

#' @noRd
validate_non_genetic <- function(non_genetic) {
  if (is.null(non_genetic)) {
    return(list())
  }
  if (inherits(non_genetic, "nongenetic_feature")) {
    non_genetic <- list(non_genetic)
  }
  if (!is.list(non_genetic)) {
    stop("`non_genetic` must be NULL, a 'nongenetic_feature', or a ",
         "list of them.", call. = FALSE)
  }
  bad <- !vapply(non_genetic, inherits, logical(1),
                 what = "nongenetic_feature")
  if (any(bad)) {
    stop("all entries of `non_genetic` must be 'nongenetic_feature' ",
         "objects; non-conforming at positions: ",
         paste(which(bad), collapse = ", "), ".", call. = FALSE)
  }
  non_genetic
}

#' @export
print.mp_case <- function(x, ...) {
  cat("<mp_case>\n")
  if (x$n_markers > 0L) {
    cat("  genetic    : ", x$n_markers, " marker(s) [",
        paste(utils::head(x$genetic$markers, 6L), collapse = ", "),
        if (x$n_markers > 6L) ", ..." else "", "]\n", sep = "")
    cat("  mutation   : ", x$genetic$mutation$model, "\n", sep = "")
  } else {
    cat("  genetic    : none\n")
  }
  if (x$n_features > 0L) {
    types <- vapply(x$non_genetic, function(f) f$type, character(1))
    cat("  non-genetic: ", x$n_features, " feature(s) [",
        paste(types, collapse = ", "), "]\n", sep = "")
  } else {
    cat("  non-genetic: none\n")
  }
  cat("  combine    : ", x$combine, "\n", sep = "")
  invisible(x)
}
