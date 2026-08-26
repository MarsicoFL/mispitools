## F4.4 — Public R API for the composed (whole-profile) LR distribution.
##
## `lr_distribution()` is the model-aware constructor promised by the
## F4.3 S3 layer: it turns a profile of `marker_model` objects into a
## single `lr_dist` object by computing each per-marker LR distribution
## (cpt_marker_joint_cpp_wrap -> cpp_per_marker_lr_dist) and composing
## them with the exact convolution kernel `cpp_lr_dist_compose` (F4.2).
## The resulting object reuses the `summary()` / `plot()` methods from
## F4.3; this file adds the `quantile()` method.

#' Composed likelihood-ratio distribution over a marker profile
#'
#' @description
#' Builds the exact distribution of the total \eqn{\log_{10}\mathrm{LR}}
#' for a profile of independent markers. Each marker's per-marker LR
#' distribution is computed from its joint H1 / H2 table, and the
#' distributions are convolved (conditional independence: the total
#' \eqn{\log_{10}\mathrm{LR}} is the sum of the per-marker
#' \eqn{\log_{10}\mathrm{LR}}s). The result is an [as_lr_dist()] object,
#' so [summary.lr_dist()], [plot.lr_dist()] and [quantile.lr_dist()]
#' apply directly.
#'
#' @param models A `marker_model` object or a `list` of them
#'   (see [marker_model()]). A single model is accepted and returns its
#'   own per-marker distribution unchanged. List names, when present and
#'   non-empty, are recorded in the `"markers"` attribute.
#' @param poi Optional character scalar naming the person of interest,
#'   applied to every model. When `NULL` (default) the POI is resolved
#'   per model by the joint engine (see [per_marker_kl()]).
#' @param method Composition method. `"exact"` (default) performs the
#'   exact sparse convolution and is bit-for-bit with the reference
#'   engine. `"grid"` projects each distribution onto a regular lattice
#'   before convolving (mass- and mean-preserving, with a shape error
#'   that shrinks with `grid_points`); it requires finite support and is
#'   useful for very long profiles.
#' @param grid_points Integer number of lattice points used when
#'   `method = "grid"` (ignored for `"exact"`). Default `512`.
#'
#' @return An object of class `lr_dist` (a data.frame with columns
#'   `log10_lr`, `p_h1`, `p_h2`), with attributes `"markers"` (the
#'   marker identifiers in input order), `"n_markers"`, `"method"` and
#'   `"poi"` (the resolved POI; a vector if it differs across models).
#'
#' @details
#' Markers are assumed conditionally independent given the hypothesis
#' (no linkage); linked-marker composition arrives with the F5 linkage
#' work. `+Inf` / `-Inf` atoms (which arise under `mutation = "none"`
#' when a hypothesis assigns zero mass to a state the other supports)
#' propagate through the exact convolution and are reported by
#' [summary.lr_dist()]; `method = "grid"` rejects infinite support.
#'
#' Two atoms carry the same `log10` LR when they agree to a relative
#' tolerance of 1e-12, rather than bit for bit. A single atom is a real
#' number that the engine can reach by more than one arithmetic route, and
#' on a platform that contracts `a * b + c` into a fused multiply-add the
#' routes differ in the last bits; grouping by exact equality would make the
#' size of the support depend on the compiler. The tolerance sits four
#' orders of magnitude above that rounding noise and four below the closest
#' genuinely distinct pair observed, so it merges duplicates and nothing
#' else.
#'
#' @seealso [marker_model()], [as_lr_dist()], [summary.lr_dist()],
#'   [plot.lr_dist()], [quantile.lr_dist()], [per_marker_kl_profile()].
#'
#' @export
#' @examples
#' if (requireNamespace("pedtools", quietly = TRUE)) {
#'   ped <- pedtools::nuclearPed(1)
#'   models <- list(
#'     D3 = marker_model(ped, "D3", c("15" = 0.4, "16" = 0.6),
#'                       mutation = list(model = "equal", rate = 1e-3)),
#'     vWA = marker_model(ped, "vWA", c("a" = 0.2, "b" = 0.3, "c" = 0.5),
#'                        mutation = list(model = "equal", rate = 1e-3))
#'   )
#'   d <- lr_distribution(models)
#'   summary(d)
#'   quantile(d)
#' }
lr_distribution <- function(models, poi = NULL,
                            method = c("exact", "grid"),
                            grid_points = 512L) {
  method <- match.arg(method)

  if (inherits(models, "marker_model")) {
    models <- list(models)
  }
  if (!is.list(models)) {
    stop("`models` must be a 'marker_model' object or a list of them.",
         call. = FALSE)
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
  if (length(grid_points) != 1L || is.na(grid_points) ||
      !is.numeric(grid_points) || grid_points < 2) {
    stop("`grid_points` must be a single integer >= 2.", call. = FALSE)
  }

  n <- length(models)
  dists <- vector("list", n)
  marker_ids <- character(n)
  pois <- character(n)
  for (i in seq_len(n)) {
    joint <- cpt_marker_joint_cpp_wrap(models[[i]], poi = poi)
    d <- cpp_per_marker_lr_dist(joint$P_H1, joint$P_H2, aggregate = TRUE)
    dists[[i]] <- list(log10_lr = d$log10_lr,
                       p_h1 = d$p_h1,
                       p_h2 = d$p_h2)
    marker_ids[i] <- models[[i]]$marker_id
    pois[i] <- attr(joint, "poi")
  }

  nm <- names(models)
  if (!is.null(nm)) {
    use_name <- nzchar(nm) & !is.na(nm)
    if (any(use_name)) marker_ids[use_name] <- nm[use_name]
  }

  comp <- cpp_lr_dist_compose(dists, method = method,
                              grid_points = as.integer(grid_points))

  poi_attr <- if (!is.null(poi)) {
    poi
  } else if (length(unique(pois)) == 1L) {
    pois[[1L]]
  } else {
    pois
  }

  out <- new_lr_dist(comp$log10_lr, comp$p_h1, comp$p_h2,
                     marker_id = NA_character_, poi = poi_attr)
  attr(out, "markers") <- marker_ids
  attr(out, "n_markers") <- n
  attr(out, "method") <- method
  out
}

#' Quantiles of a likelihood-ratio distribution
#'
#' @description
#' Discrete inverse-CDF quantiles of \eqn{\log_{10}\mathrm{LR}} under
#' one of the two hypotheses, using R's `quantile` type-1 definition
#' \eqn{Q(p) = \inf\{x : F(x) \ge p\}} over the active support.
#'
#' @param x An `lr_dist` object (see [as_lr_dist()] / [lr_distribution()]).
#' @param probs Numeric vector of probabilities in \eqn{[0, 1]}.
#'   Defaults to the quartiles `c(0, 0.25, 0.5, 0.75, 1)`.
#' @param under_h1 Logical. If `TRUE` (default) quantiles are taken
#'   under \eqn{H_1} (weights `p_h1`); if `FALSE`, under \eqn{H_2}
#'   (weights `p_h2`).
#' @param ... Unused.
#'
#' @return A named numeric vector of quantiles, one per `probs` entry.
#'
#' @details
#' Atoms with zero mass under the requested hypothesis are dropped
#' before the cumulative distribution is formed. A `+Inf` (or `-Inf`)
#' atom is a legitimate quantile value when the requested probability
#' falls in its mass. An error is raised if no atom has positive mass
#' under the requested hypothesis.
#'
#' @seealso [lr_distribution()], [summary.lr_dist()], [plot.lr_dist()]
#'
#' @importFrom stats quantile
#' @export
#' @examples
#' d <- as_lr_dist(data.frame(
#'   log10_lr = c(-1, 0, 2),
#'   p_h1 = c(0.1, 0.3, 0.6),
#'   p_h2 = c(0.6, 0.3, 0.1)
#' ))
#' quantile(d)
#' quantile(d, probs = c(0.05, 0.95), under_h1 = FALSE)
quantile.lr_dist <- function(x, probs = c(0, 0.25, 0.5, 0.75, 1),
                             under_h1 = TRUE, ...) {
  if (!is.numeric(probs) || length(probs) == 0L || anyNA(probs) ||
      any(probs < 0) || any(probs > 1)) {
    stop("`probs` must be a numeric vector with values in [0, 1].",
         call. = FALSE)
  }
  if (length(under_h1) != 1L || is.na(under_h1) ||
      !is.logical(under_h1)) {
    stop("`under_h1` must be a single logical value.", call. = FALSE)
  }
  q <- cpp_lr_dist_quantile(x$log10_lr, x$p_h1, x$p_h2,
                            as.numeric(probs), under_h1 = under_h1)
  names(q) <- paste0(format(100 * probs, trim = TRUE), "%")
  q
}
