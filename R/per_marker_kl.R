## F3.2 — Public R API for per-marker bidirectional KL.
##
## Thin wrapper layer over `cpt_marker_joint_cpp_wrap()` + `cpp_per_marker_kl()`.
## Mirrors the column set of the F1.6 reference engine `per_marker_kl_R()` and
## extends it with the KLde-style absolute-continuity diagnostics surfaced by
## the C++ kernel (F3.1). As of F3.4c the profile entry point routes through
## the single-trip core batch path (F3.4a primitive + F3.4b mutation-matrix
## cache); the R-level per-marker loop survives only as a fallback for
## heterogeneous topology, linkage, or unsupported mutation models.

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
#' order. When the profile shares one pedigree topology, uses no linkage,
#' and every mutation model is wired to the C++ backend
#' (`none` / `equal` / `stepwise`), evaluation routes through a single
#' cross-marker C++ batch call with a shared mutation-matrix cache;
#' otherwise it falls back to a per-marker R-level loop.
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

  ## F3.4 — batch path when all models share pedigree topology, none use
  ## linkage, and the mutation kinds are all wired to the C++ backend
  ## (none / equal / stepwise). Falls through to the per-marker loop
  ## otherwise.
  batch <- per_marker_kl_profile_batch(models, poi)
  out <- if (!is.null(batch)) {
    batch
  } else {
    rows <- vector("list", length(models))
    pois <- character(length(models))
    for (i in seq_along(models)) {
      r <- per_marker_kl(models[[i]], poi = poi)
      pois[i] <- attr(r, "poi")
      rows[[i]] <- r
    }
    out0 <- do.call(rbind, c(rows, list(make.row.names = FALSE)))
    attr(out0, "poi") <- if (!is.null(poi)) poi else pois
    out0
  }

  nm <- names(models)
  if (!is.null(nm)) {
    use_name <- nzchar(nm) & !is.na(nm)
    if (any(use_name)) {
      out$marker[use_name] <- nm[use_name]
    }
  }

  out
}

#' @noRd
per_marker_kl_profile_batch <- function(models, poi) {
  ## Returns a data.frame in the per_marker_kl()-output shape when the
  ## batch path applies, or NULL when the caller must use the R-level
  ## fallback (heterogeneous topology, linkage, unsupported mutation).
  m1 <- models[[1L]]
  if (!is.null(m1$linkage)) return(NULL)
  if (!m1$mutation$model %in% c("none", "equal", "stepwise")) return(NULL)
  if (length(models) > 1L) {
    for (i in 2:length(models)) {
      m <- models[[i]]
      if (!is.null(m$linkage)) return(NULL)
      if (!m$mutation$model %in% c("none", "equal", "stepwise")) return(NULL)
      if (!identical(m$ped, m1$ped)) return(NULL)
    }
  }

  if (!requireNamespace("pedtools", quietly = TRUE)) {
    stop("Package 'pedtools' is required.", call. = FALSE)
  }

  ped <- m1$ped
  members <- as.character(labels(ped))
  n <- length(members)
  name_to_idx <- stats::setNames(seq_len(n) - 1L, members)
  father <- integer(n)
  mother <- integer(n)
  for (i in seq_len(n)) {
    fa <- as.character(pedtools::father(ped, id = members[i]))
    mo <- as.character(pedtools::mother(ped, id = members[i]))
    father[i] <- if (length(fa) == 0L || !nzchar(fa) || is.na(fa)) -1L else
      unname(name_to_idx[[fa]])
    mother[i] <- if (length(mo) == 0L || !nzchar(mo) || is.na(mo)) -1L else
      unname(name_to_idx[[mo]])
  }
  poi_id <- resolve_poi(ped, poi)
  poi_idx <- unname(name_to_idx[[poi_id]])

  N <- length(models)
  freqs_list <- vector("list", N)
  labels_list <- vector("list", N)
  kind_vec <- integer(N)
  rate_vec <- numeric(N)
  range_vec <- numeric(N)
  marker_ids <- character(N)

  for (i in seq_len(N)) {
    m <- models[[i]]
    freqs_list[[i]] <- as.numeric(unname(m$freqs))
    K <- length(m$alleles)
    kind_vec[i] <- switch(m$mutation$model, none = 0L, equal = 1L, stepwise = 2L)
    rate_vec[i] <- if (m$mutation$model == "none") 0.0 else as.numeric(m$mutation$rate)
    range_vec[i] <- if (m$mutation$model == "stepwise")
      as.numeric(m$mutation$ratio) else 0.0
    if (m$mutation$model == "stepwise") {
      s <- suppressWarnings(as.numeric(m$alleles))
      if (anyNA(s)) {
        stop("Stepwise mutation requires numeric allele labels; got non-",
             "numeric: ", paste(m$alleles[is.na(s)], collapse = ", "),
             call. = FALSE)
      }
      labels_list[[i]] <- s
    } else {
      labels_list[[i]] <- rep(NA_real_, K)
    }
    marker_ids[i] <- m$marker_id
  }

  res <- cpp_per_marker_kl_batch(
    father = as.integer(father),
    mother = as.integer(mother),
    poi = as.integer(poi_idx),
    freqs_list = freqs_list,
    mutation_kind = kind_vec,
    mutation_rate = rate_vec,
    mutation_range = range_vec,
    numeric_labels_list = labels_list
  )

  out <- data.frame(
    marker = marker_ids,
    e_log10_lr_h1 = res$e_log10_lr_h1,
    e_log10_lr_h2 = res$e_log10_lr_h2,
    kl_h1_to_h2  = res$kl_h1_to_h2,
    kl_h2_to_h1  = res$kl_h2_to_h1,
    abs_cont_violations_h1 = as.integer(res$abs_cont_violations_h1),
    abs_cont_violations_h2 = as.integer(res$abs_cont_violations_h2),
    mass_violations_h1 = res$mass_violations_h1,
    mass_violations_h2 = res$mass_violations_h2,
    stringsAsFactors = FALSE
  )
  attr(out, "poi") <- if (!is.null(poi)) poi else rep(poi_id, N)
  attr(out, "cache_hits")   <- as.integer(res$cache_hits)
  attr(out, "cache_misses") <- as.integer(res$cache_misses)
  out
}
