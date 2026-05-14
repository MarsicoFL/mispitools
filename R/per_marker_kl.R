## F3.2 — Public R API for per-marker bidirectional KL.
##
## Thin wrapper layer over `cpt_marker_joint_cpp_wrap()` + `cpp_per_marker_kl()`.
## Mirrors the column set of the F1.6 reference engine `per_marker_kl_R()` and
## extends it with the KLde-style absolute-continuity diagnostics surfaced by
## the C++ kernel (F3.1). Core-side vectorisation across markers + mutation
## matrix caching lands in F3.4; here the profile entry point loops in R.

#' Per-marker bidirectional Kullback-Leibler divergence and expected log10 LR
#'
#' @description
#' Computes, for a single marker model, the expected `log10` likelihood ratio
#' under both hypotheses and the bidirectional Kullback-Leibler divergence
#' between the joint H1 and H2 distributions over pedigree typings. Drives the
#' C++ engine introduced in F2 / F3.1; returns the same numerical content as
#' the internal reference engine to bit-equivalent precision.
#'
#' Boundary convention follows the reference engine:
#' \itemize{
#'   \item `P_H1(g) > 0` and `P_H2(g) = 0` contribute `+Inf` to
#'         `kl_h1_to_h2` and `e_log10_lr_h1`.
#'   \item `P_H1(g) = 0` and `P_H2(g) > 0` contribute `+Inf` to
#'         `kl_h2_to_h1` and `-Inf` to `e_log10_lr_h2`.
#'   \item States with `P_H1 = P_H2 = 0` are skipped (already filtered by the
#'         joint engine).
#' }
#' Under `mutation = "none"` the H2 distribution typically spans
#' Mendelian-incompatible states with zero mass under H1, so
#' `kl_h2_to_h1 = +Inf`. The `abs_cont_violations_*` and
#' `mass_violations_*` columns quantify how much probability mass falls on
#' violating states in each direction.
#'
#' @param model A `marker_model` object (see [marker_model()]).
#' @param poi Optional character scalar naming the person of interest in the
#'   pedigree. When `NULL` (default), the engine resolves the POI by picking
#'   the last untyped non-founder, falling back to the last non-founder, and
#'   finally to the last member.
#'
#' @return A one-row `data.frame` with columns
#'   \describe{
#'     \item{marker}{Marker identifier (`character`).}
#'     \item{e_log10_lr_h1}{Expected `log10(LR)` under H1.}
#'     \item{e_log10_lr_h2}{Expected `log10(LR)` under H2.}
#'     \item{kl_h1_to_h2}{KL(H1 || H2) in nats.}
#'     \item{kl_h2_to_h1}{KL(H2 || H1) in nats.}
#'     \item{abs_cont_violations_h1}{Number of joint states with positive H2
#'           mass but zero H1 mass (integer).}
#'     \item{abs_cont_violations_h2}{Number of joint states with positive H1
#'           mass but zero H2 mass (integer).}
#'     \item{mass_violations_h1}{Total H2 mass on states violating H1 absolute
#'           continuity (numeric).}
#'     \item{mass_violations_h2}{Total H1 mass on states violating H2 absolute
#'           continuity (numeric).}
#'   }
#' Attribute `"poi"` records the resolved POI identifier.
#'
#' @seealso [marker_model()], [per_marker_kl_profile()].
#'
#' @export
#' @examples
#' if (requireNamespace("pedtools", quietly = TRUE)) {
#'   ped <- pedtools::nuclearPed(1)
#'   freqs <- c("a" = 0.3, "b" = 0.5, "c" = 0.2)
#'   mm <- marker_model(ped, "M1", freqs,
#'                      mutation = list(model = "equal", rate = 0.005))
#'   per_marker_kl(mm)
#' }
per_marker_kl <- function(model, poi = NULL) {
  if (!inherits(model, "marker_model")) {
    stop("`model` must be a 'marker_model' object.", call. = FALSE)
  }
  joint <- cpt_marker_joint_cpp_wrap(model, poi = poi)
  kl <- cpp_per_marker_kl(joint$P_H1, joint$P_H2)

  out <- data.frame(
    marker = model$marker_id,
    e_log10_lr_h1 = kl$e_log10_lr_h1,
    e_log10_lr_h2 = kl$e_log10_lr_h2,
    kl_h1_to_h2 = kl$kl_h1_to_h2,
    kl_h2_to_h1 = kl$kl_h2_to_h1,
    abs_cont_violations_h1 = as.integer(kl$abs_cont_violations_h1),
    abs_cont_violations_h2 = as.integer(kl$abs_cont_violations_h2),
    mass_violations_h1 = kl$mass_violations_h1,
    mass_violations_h2 = kl$mass_violations_h2,
    stringsAsFactors = FALSE
  )
  attr(out, "poi") <- attr(joint, "poi")
  out
}

#' Per-marker bidirectional KL across a marker profile
#'
#' @description
#' Vectorised wrapper around [per_marker_kl()] for a list of `marker_model`
#' objects. Returns a `data.frame` with one row per input model in input
#' order. F3.2 loops at the R level; the cross-marker C++ batch entry point
#' (with mutation matrix caching) arrives in F3.4.
#'
#' All models must share the same pedigree topology if `poi` is supplied as
#' a scalar; otherwise the POI is resolved independently for each model.
#'
#' @param models A `list` of `marker_model` objects. Names of the list, when
#'   present and non-empty, override the per-model `marker_id` in the output
#'   `marker` column.
#' @param poi Optional character scalar applied to every model, or `NULL`
#'   (default) for per-model resolution. Pass a vector by mapping the loop
#'   yourself if you need heterogeneous POIs.
#'
#' @return A `data.frame` with the same columns as [per_marker_kl()] and
#'   `nrow(out) == length(models)`. When `poi` is a scalar, attribute
#'   `"poi"` carries that value; otherwise the resolved POIs appear as
#'   attribute `"poi"` of length `length(models)`.
#'
#' @seealso [per_marker_kl()].
#'
#' @export
#' @examples
#' if (requireNamespace("pedtools", quietly = TRUE)) {
#'   ped <- pedtools::nuclearPed(1)
#'   models <- list(
#'     M1 = marker_model(ped, "M1", c("a" = 0.4, "b" = 0.6),
#'                       mutation = list(model = "equal", rate = 1e-3)),
#'     M2 = marker_model(ped, "M2", c("a" = 0.2, "b" = 0.3, "c" = 0.5),
#'                       mutation = list(model = "equal", rate = 1e-3))
#'   )
#'   per_marker_kl_profile(models)
#' }
per_marker_kl_profile <- function(models, poi = NULL) {
  if (!is.list(models)) {
    stop("`models` must be a list of 'marker_model' objects.", call. = FALSE)
  }
  if (length(models) == 0L) {
    stop("`models` must contain at least one 'marker_model' object.",
         call. = FALSE)
  }
  bad <- !vapply(models, inherits, logical(1), what = "marker_model")
  if (any(bad)) {
    stop(sprintf(
      "All entries of `models` must be 'marker_model' objects; non-conforming entries at positions: %s.",
      paste(which(bad), collapse = ", ")), call. = FALSE)
  }
  if (!is.null(poi)) {
    if (!is.character(poi) || length(poi) != 1L || is.na(poi) ||
        !nzchar(poi)) {
      stop("`poi` must be NULL or a single non-empty character string.",
           call. = FALSE)
    }
  }

  rows <- vector("list", length(models))
  pois <- character(length(models))
  for (i in seq_along(models)) {
    r <- per_marker_kl(models[[i]], poi = poi)
    pois[i] <- attr(r, "poi")
    rows[[i]] <- r
  }
  out <- do.call(rbind, c(rows, list(make.row.names = FALSE)))

  nm <- names(models)
  if (!is.null(nm)) {
    use_name <- nzchar(nm) & !is.na(nm)
    if (any(use_name)) {
      out$marker[use_name] <- nm[use_name]
    }
  }

  attr(out, "poi") <- if (!is.null(poi)) poi else pois
  out
}
