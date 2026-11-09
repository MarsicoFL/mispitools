
#' Evaluate all evidence in a missing-person case
#'
#' @description
#' Runs the unified evaluation pipeline on an [mp_case()] object and
#' returns every reportable quantity in a single object: the
#' per-feature Kullback-Leibler / expected-`log10`-LR table (genetic
#' and non-genetic in one frame), the per-feature LR distributions,
#' the combined LR distribution, the expected combined `log10` LR
#' under each hypothesis, the decision-theoretic quantities
#' (threshold, FPR, FNR), and the fragility diagnostics (concentration
#' index and leave-one-out).
#'
#' @param case An [mp_case()] object.
#' @param poi Optional character scalar naming the person of interest
#'   in the pedigree, passed to the genetic engine. `NULL` (default)
#'   lets the joint engine resolve the POI.
#' @param n_threads Optional integer thread count for the engine.
#'   `NULL` (default) resolves the CRAN-safe policy via the
#'   `mispitools.n_threads` option (see the parallel build notes). The
#'   resolved value is recorded on the result.
#'
#' @return An object of class `"evidence_result"`: a list with
#'   \describe{
#'     \item{`per_feature_kl`}{One row per feature with columns
#'       `feature`, `kind` (`"genetic"` / `"non_genetic"`),
#'       `e_log10_lr_h1`, `e_log10_lr_h2`, `kl_h1_to_h2`,
#'       `kl_h2_to_h1`.}
#'     \item{`per_feature_lr_dist`}{Named list of per-feature sparse
#'       `(log10_lr, p_h1, p_h2)` data.frames.}
#'     \item{`combined_lr_dist`}{The combined distribution as an
#'       [as_lr_dist()] object.}
#'     \item{`expected_log10_lr_h1`, `expected_log10_lr_h2`}{Means of
#'       the combined `log10` LR under each hypothesis.}
#'     \item{`decision_quantities`}{List with `threshold`, `fpr`,
#'       `fnr`, `tpr`, `tnr`, `mcc`, `auc`.}
#'     \item{`fragility`}{List with `concentration`, `leave_one_out`
#'       and (when all weights are finite) a `report`.}
#'     \item{`combine`, `n_threads`}{The combination mode used and the
#'       resolved thread count.}
#'   }
#'
#' @details
#' Genetic per-marker KL routes through [per_marker_kl_profile()]
#' (single-trip C++ batch when the profile is homogeneous);
#' non-genetic features go through the F6 kernel. The combined
#' distribution uses `cpp_evidence_combine()` with the case's
#' `combine` mode, so under `"independent"` it is bit-for-bit the
#' exact convolution and under `"markov_se"` it follows the
#' Egeland-Marsico (2026) chain. Under `mutation = "none"` per-feature
#' expectations can be infinite; the concentration index then uses the
#' finite weights and the `report` is omitted.
#'
#' @seealso [mp_case()], [per_marker_kl_profile()],
#'   [lr_distribution()], [fragility_report()], [as_lr_dist()].
#'
#' @export
#' @examples
#' if (requireNamespace("pedtools", quietly = TRUE)) {
#'   ped <- pedtools::nuclearPed(1)
#'   case <- mp_case(
#'     genetic = list(
#'       ped = ped,
#'       freqs = Argentina,
#'       markers = c("D3S1358", "THO1", "VWA"),
#'       mutation = list(model = "equal", rate = 1e-3)
#'     ),
#'     non_genetic = list(
#'       sex = nongenetic_feature(
#'         type = "sex", observed = "F",
#'         db_or_freqs = c(F = 0.5, M = 0.5), error = 0.05
#'       )
#'     )
#'   )
#'   res <- evaluate_evidence(case)
#'   res$per_feature_kl
#'   res$expected_log10_lr_h1
#' }
evaluate_evidence <- function(case, poi = NULL, n_threads = NULL) {
  if (!inherits(case, "mp_case")) {
    stop("`case` must be an 'mp_case' object (see mp_case()).",
         call. = FALSE)
  }
  if (!is.null(poi)) {
    if (!is.character(poi) || length(poi) != 1L || is.na(poi) ||
        !nzchar(poi)) {
      stop("`poi` must be NULL or a single non-empty character string.",
           call. = FALSE)
    }
  }

  resolved_threads <- mispi_n_threads(n_threads)
  old_opt <- options(mispitools.n_threads = resolved_threads)
  on.exit(options(old_opt), add = TRUE)

  models <- if (is.null(case$genetic)) list() else case$genetic$models
  features <- case$non_genetic

  kl_rows <- list()
  dists <- list()
  dist_names <- character(0)

  if (length(models) > 0L) {
    kl_g <- per_marker_kl_profile(models, poi = poi)
    kl_rows[["genetic"]] <- data.frame(
      feature = kl_g$marker,
      kind = "genetic",
      e_log10_lr_h1 = kl_g$e_log10_lr_h1,
      e_log10_lr_h2 = kl_g$e_log10_lr_h2,
      kl_h1_to_h2 = kl_g$kl_h1_to_h2,
      kl_h2_to_h1 = kl_g$kl_h2_to_h1,
      stringsAsFactors = FALSE
    )
    for (i in seq_along(models)) {
      joint <- cpt_marker_joint_cpp_wrap(models[[i]], poi = poi)
      d <- cpp_per_marker_lr_dist(joint$P_H1, joint$P_H2,
                                  aggregate = TRUE)
      nm <- names(models)[i] %||% models[[i]]$marker_id
      dists[[length(dists) + 1L]] <- list(log10_lr = d$log10_lr,
                                          p_h1 = d$p_h1,
                                          p_h2 = d$p_h2)
      dist_names <- c(dist_names, nm)
    }
  }

  if (length(features) > 0L) {
    fnames <- names(features)
    ng_kl <- vector("list", length(features))
    for (i in seq_along(features)) {
      fk <- per_feature_kl_cpp_wrap(features[[i]])
      lab <- if (!is.null(fnames) && nzchar(fnames[i])) fnames[i] else
        fk$feature
      ng_kl[[i]] <- data.frame(
        feature = lab, kind = "non_genetic",
        e_log10_lr_h1 = fk$e_log10_lr_h1,
        e_log10_lr_h2 = fk$e_log10_lr_h2,
        kl_h1_to_h2 = fk$kl_h1_to_h2,
        kl_h2_to_h1 = fk$kl_h2_to_h1,
        stringsAsFactors = FALSE
      )
      fd <- per_feature_lr_dist_cpp_wrap(features[[i]], aggregate = TRUE)
      dists[[length(dists) + 1L]] <- list(log10_lr = fd$log10_lr,
                                          p_h1 = fd$p_h1,
                                          p_h2 = fd$p_h2)
      dist_names <- c(dist_names, lab)
    }
    kl_rows[["non_genetic"]] <-
      do.call(rbind, c(ng_kl, list(make.row.names = FALSE)))
  }

  per_feature_kl <- do.call(rbind,
                            c(kl_rows, list(make.row.names = FALSE)))

  comp <- cpp_evidence_combine(dists, mode = case$combine,
                               transition = case$transition)
  combined <- new_lr_dist(comp$log10_lr, comp$p_h1, comp$p_h2,
                          marker_id = NA_character_, poi = poi)
  attr(combined, "features") <- dist_names
  attr(combined, "n_features") <- length(dists)

  names(dists) <- dist_names
  per_feature_lr_dist <- lapply(dists, function(d) {
    data.frame(log10_lr = d$log10_lr, p_h1 = d$p_h1, p_h2 = d$p_h2,
               stringsAsFactors = FALSE)
  })

  s <- cpp_lr_dist_summary(combined$log10_lr, combined$p_h1,
                           combined$p_h2)
  ct <- cpp_lr_dist_choose_threshold(combined$log10_lr, combined$p_h1,
                                     combined$p_h2, weight = 10.0)
  rates <- cpp_lr_dist_decision_rates(combined$log10_lr, combined$p_h1,
                                      combined$p_h2,
                                      threshold = ct$threshold)
  roc <- cpp_lr_dist_roc(combined$log10_lr, combined$p_h1,
                         combined$p_h2)
  decision_quantities <- list(
    threshold = ct$threshold,
    fpr = rates$fpr,
    fnr = rates$fnr,
    tpr = rates$tpr,
    tnr = rates$tnr,
    mcc = rates$mcc,
    auc = roc$auc
  )

  w <- stats::setNames(per_feature_kl$e_log10_lr_h1,
                       per_feature_kl$feature)
  finite_w <- w[is.finite(w)]
  conc <- concentration_index_positive(unname(finite_w))
  if (all(is.finite(w))) {
    lrs <- stats::setNames(10^w, names(w))
    loo <- leave_one_out(lrs)
    rep_obj <- fragility_report(lrs)
  } else {
    loo <- data.frame(
      feature = names(w),
      e_log10_lr_h1 = unname(w),
      total_without = sum(finite_w) - ifelse(is.finite(w), w, 0),
      stringsAsFactors = FALSE
    )
    rep_obj <- NULL
  }
  fragility <- list(concentration = conc,
                    leave_one_out = loo,
                    report = rep_obj)

  structure(
    list(
      per_feature_kl = per_feature_kl,
      per_feature_lr_dist = per_feature_lr_dist,
      combined_lr_dist = combined,
      expected_log10_lr_h1 = s$mean_h1,
      expected_log10_lr_h2 = s$mean_h2,
      decision_quantities = decision_quantities,
      fragility = fragility,
      combine = case$combine,
      n_threads = resolved_threads
    ),
    class = "evidence_result"
  )
}

#' @export
print.evidence_result <- function(x, ...) {
  cat("<evidence_result>\n")
  cat("  features            : ", nrow(x$per_feature_kl),
      " (combine: ", x$combine, ")\n", sep = "")
  cat("  E[log10 LR | H1]    : ",
      format(x$expected_log10_lr_h1, digits = 5), "\n", sep = "")
  cat("  E[log10 LR | H2]    : ",
      format(x$expected_log10_lr_h2, digits = 5), "\n", sep = "")
  cat("  decision threshold  : ",
      format(x$decision_quantities$threshold, digits = 5),
      "  (FPR ", format(x$decision_quantities$fpr, digits = 3),
      ", FNR ", format(x$decision_quantities$fnr, digits = 3),
      ")\n", sep = "")
  cat("  concentration C_W+  : ",
      format(x$fragility$concentration, digits = 4), "\n", sep = "")
  invisible(x)
}

#' @export
summary.evidence_result <- function(object, ...) {
  cat("Evidence evaluation summary\n")
  cat("Per-feature expected log10 LR (under H1):\n")
  tab <- object$per_feature_kl[, c("feature", "kind",
                                   "e_log10_lr_h1", "e_log10_lr_h2")]
  print(tab, row.names = FALSE)
  cat("\nCombined: E[log10 LR | H1] = ",
      format(object$expected_log10_lr_h1, digits = 5),
      ", E[log10 LR | H2] = ",
      format(object$expected_log10_lr_h2, digits = 5), "\n", sep = "")
  cat("Decision: threshold = ",
      format(object$decision_quantities$threshold, digits = 5),
      ", AUC = ", format(object$decision_quantities$auc, digits = 5),
      "\n", sep = "")
  cat("Fragility: C_W+ = ",
      format(object$fragility$concentration, digits = 4), "\n",
      sep = "")
  invisible(object)
}
