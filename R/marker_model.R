#' Marker model (S3)
#'
#' @description
#' Constructor for a single-marker forensic genetic model. Bundles the
#' pedigree, the marker identifier, the population allele frequencies, the
#' mutation model, and (optionally) a linkage descriptor. The resulting
#' object is the unit consumed by the per-marker engines (`per_marker_kl`,
#' `per_marker_lr_dist`) that arrive in later milestones.
#'
#' This is the public entry point for the F1 reference engine. The
#' constructor only validates and stores its inputs; the actual joint CPT
#' construction is done downstream.
#'
#' @param ped A `pedtools::ped` object describing the pedigree. Singletons
#'   and `pedList` objects are not accepted at this stage.
#' @param marker_id Character scalar. Identifier used to label the marker
#'   (e.g. `"D3S1358"`). When `ped` already carries a marker with this
#'   name, downstream engines will read the genotypes of typed individuals
#'   from it; otherwise the model represents the prospective situation
#'   prior to typing.
#' @param freqs Named numeric vector of population allele frequencies. Names
#'   are allele labels, values must lie in `[0, 1]` and sum to 1 within
#'   `tol`. At least two alleles are required.
#' @param mutation List describing the mutation model. Required field
#'   `model` is one of `"none"`, `"equal"`, `"stepwise"`, `"asymmetric"`.
#'   For `"none"` no further fields are needed; for the other models a
#'   numeric `rate` in `[0, 1)` is required. `"stepwise"` additionally
#'   accepts `ratio` (geometric step ratio in `(0, 1)`); `"asymmetric"`
#'   additionally accepts `ratio` and `bias` (`u` parameter, in `[0, 1]`).
#'   Defaults to `list(model = "none", rate = 0)`.
#' @param linkage Either `NULL` (default; unlinked marker) or a list with
#'   fields `partner` (character, the marker identifier of the linked
#'   partner) and `theta` (recombination fraction in `[0, 0.5]`). Linkage
#'   support is wired in milestone F5; F1 only validates the structure.
#' @param tol Numeric tolerance used when checking that `freqs` sums to 1.
#'   Defaults to `1e-6`.
#'
#' @return An object of class `"marker_model"` carrying the validated
#'   inputs as components `ped`, `marker_id`, `freqs`, `mutation`,
#'   `linkage`, and `alleles` (`names(freqs)`).
#'
#' @export
#' @examples
#' if (requireNamespace("pedtools", quietly = TRUE)) {
#'   ped <- pedtools::nuclearPed(1)
#'   freqs <- c("12" = 0.2, "13" = 0.3, "14" = 0.5)
#'   mm <- marker_model(
#'     ped = ped,
#'     marker_id = "M1",
#'     freqs = freqs,
#'     mutation = list(model = "equal", rate = 1e-3)
#'   )
#'   print(mm)
#' }
marker_model <- function(ped,
                         marker_id,
                         freqs,
                         mutation = list(model = "none", rate = 0),
                         linkage = NULL,
                         tol = 1e-6) {
  validate_ped(ped)
  marker_id <- validate_marker_id(marker_id)
  freqs <- validate_freqs(freqs, tol = tol)
  mutation <- validate_mutation(mutation)
  linkage <- validate_linkage(linkage)

  structure(
    list(
      ped = ped,
      marker_id = marker_id,
      freqs = freqs,
      alleles = names(freqs),
      mutation = mutation,
      linkage = linkage
    ),
    class = "marker_model"
  )
}

#' @noRd
validate_ped <- function(ped) {
  if (!requireNamespace("pedtools", quietly = TRUE)) {
    stop("Package 'pedtools' is required for marker_model().", call. = FALSE)
  }
  if (!pedtools::is.ped(ped)) {
    stop("`ped` must be a 'ped' object from package pedtools (singletons ",
         "and pedLists are not accepted in this constructor).",
         call. = FALSE)
  }
  invisible(ped)
}

#' @noRd
validate_marker_id <- function(marker_id) {
  if (!is.character(marker_id) || length(marker_id) != 1L || is.na(marker_id)) {
    stop("`marker_id` must be a single non-NA character string.",
         call. = FALSE)
  }
  if (!nzchar(marker_id)) {
    stop("`marker_id` must be a non-empty string.", call. = FALSE)
  }
  marker_id
}

#' @noRd
validate_freqs <- function(freqs, tol = 1e-6) {
  if (!is.numeric(freqs) || is.array(freqs)) {
    stop("`freqs` must be a numeric vector of allele frequencies.",
         call. = FALSE)
  }
  nm <- names(freqs)
  if (is.null(nm) || any(!nzchar(nm)) || any(is.na(nm))) {
    stop("`freqs` must be named (one name per allele).", call. = FALSE)
  }
  if (anyDuplicated(nm)) {
    stop("`freqs` must not contain duplicated allele labels.",
         call. = FALSE)
  }
  if (length(freqs) < 2L) {
    stop("`freqs` must contain at least two alleles.", call. = FALSE)
  }
  if (anyNA(freqs)) {
    stop("`freqs` must not contain NA values.", call. = FALSE)
  }
  if (any(freqs < 0) || any(freqs > 1)) {
    stop("`freqs` entries must lie in [0, 1].", call. = FALSE)
  }
  s <- sum(freqs)
  if (abs(s - 1) > tol) {
    stop(sprintf("`freqs` must sum to 1 (within tol=%g); got sum=%.10f.",
                 tol, s), call. = FALSE)
  }
  freqs
}

#' @noRd
validate_mutation <- function(mutation) {
  if (!is.list(mutation) || is.null(names(mutation))) {
    stop("`mutation` must be a named list with at least a `model` field.",
         call. = FALSE)
  }
  if (is.null(mutation$model)) {
    stop("`mutation` must contain a `model` field.", call. = FALSE)
  }
  allowed <- c("none", "equal", "stepwise", "asymmetric")
  if (!is.character(mutation$model) || length(mutation$model) != 1L ||
      !mutation$model %in% allowed) {
    stop("`mutation$model` must be one of: ",
         paste(shQuote(allowed), collapse = ", "), ".",
         call. = FALSE)
  }

  if (mutation$model == "none") {
    if (is.null(mutation$rate)) {
      mutation$rate <- 0
    } else if (!is_scalar_finite(mutation$rate) || mutation$rate != 0) {
      stop("`mutation$rate` must be 0 when model = \"none\".",
           call. = FALSE)
    }
    return(mutation)
  }

  if (is.null(mutation$rate)) {
    stop(sprintf("`mutation$rate` is required when model = \"%s\".",
                 mutation$model), call. = FALSE)
  }
  if (!is_scalar_finite(mutation$rate) || mutation$rate < 0 || mutation$rate >= 1) {
    stop("`mutation$rate` must be a finite numeric scalar in [0, 1).",
         call. = FALSE)
  }

  if (mutation$model == "stepwise") {
    if (!is.null(mutation$ratio)) {
      if (!is_scalar_finite(mutation$ratio) || mutation$ratio <= 0 ||
          mutation$ratio >= 1) {
        stop("`mutation$ratio` must be a finite numeric scalar in (0, 1).",
             call. = FALSE)
      }
    }
  }

  if (mutation$model == "asymmetric") {
    if (!is.null(mutation$ratio)) {
      if (!is_scalar_finite(mutation$ratio) || mutation$ratio <= 0 ||
          mutation$ratio >= 1) {
        stop("`mutation$ratio` must be a finite numeric scalar in (0, 1).",
             call. = FALSE)
      }
    }
    if (!is.null(mutation$bias)) {
      if (!is_scalar_finite(mutation$bias) || mutation$bias < 0 ||
          mutation$bias > 1) {
        stop("`mutation$bias` must be a finite numeric scalar in [0, 1].",
             call. = FALSE)
      }
    }
  }

  mutation
}

#' @noRd
validate_linkage <- function(linkage) {
  if (is.null(linkage)) {
    return(NULL)
  }
  if (!is.list(linkage) || is.null(names(linkage))) {
    stop("`linkage` must be NULL or a named list.", call. = FALSE)
  }
  if (is.null(linkage$partner) || !is.character(linkage$partner) ||
      length(linkage$partner) != 1L || is.na(linkage$partner) ||
      !nzchar(linkage$partner)) {
    stop("`linkage$partner` must be a single non-empty character string.",
         call. = FALSE)
  }
  if (is.null(linkage$theta) || !is_scalar_finite(linkage$theta) ||
      linkage$theta < 0 || linkage$theta > 0.5) {
    stop("`linkage$theta` must be a finite numeric scalar in [0, 0.5].",
         call. = FALSE)
  }
  linkage
}

#' @noRd
is_scalar_finite <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x)
}

#' Print method for marker_model
#'
#' @param x A `marker_model` object.
#' @param ... Unused.
#' @return Invisibly returns `x`.
#' @export
print.marker_model <- function(x, ...) {
  cat("<marker_model>\n")
  cat("  marker_id : ", x$marker_id, "\n", sep = "")
  cat("  alleles   : ", length(x$alleles),
      " (", paste(utils::head(x$alleles, 6), collapse = ", "),
      if (length(x$alleles) > 6L) ", ..." else "",
      ")\n", sep = "")
  mut <- x$mutation
  mut_desc <- mut$model
  if (mut$model != "none") {
    mut_desc <- paste0(mut_desc, " (rate=", format(mut$rate, digits = 3))
    if (!is.null(mut$ratio)) {
      mut_desc <- paste0(mut_desc, ", ratio=", format(mut$ratio, digits = 3))
    }
    if (!is.null(mut$bias)) {
      mut_desc <- paste0(mut_desc, ", bias=", format(mut$bias, digits = 3))
    }
    mut_desc <- paste0(mut_desc, ")")
  }
  cat("  mutation  : ", mut_desc, "\n", sep = "")
  if (is.null(x$linkage)) {
    cat("  linkage   : none\n")
  } else {
    cat("  linkage   : partner=", x$linkage$partner,
        ", theta=", format(x$linkage$theta, digits = 3), "\n", sep = "")
  }
  cat("  pedigree  : ", pedtools::pedsize(x$ped), " individuals\n", sep = "")
  invisible(x)
}
