#' Calibrate pedigree-specific concentration cutoff
#'
#' Simulates the distribution of the inclusion-fragility concentration
#' index \eqn{C_W^+} under the prosecution hypothesis \eqn{H_p} for a
#' given pedigree and marker panel, and returns the chosen upper
#' quantile as a case-flagging cutoff. Profiles whose observed
#' \eqn{C_W^+} exceeds the cutoff should trigger a leave-one-out
#' review before reporting.
#'
#' @param reference A \code{pedtools::ped} object with markers attached
#'   and founder profiles simulated, as expected by
#'   \code{\link{sim_lr_genetic}}.
#' @param missing Character or numeric. ID of the missing person in the
#'   pedigree.
#' @param numsims Integer. Number of \eqn{H_p} simulations used to
#'   build the empirical \eqn{C_W^+} distribution. Default 1500.
#' @param probs Numeric in \eqn{(0, 1)}. Quantile of the \eqn{H_p}
#'   distribution to use as the cutoff. Default 0.90 (matches the
#'   paper convention).
#' @param seed Integer. Random seed passed through to
#'   \code{\link{sim_lr_genetic}}.
#' @param numCores Integer. Cores for the underlying
#'   \code{forrel::profileSim} call.
#'
#' @return A list with elements
#'   \describe{
#'     \item{\code{cutoff}}{The \code{probs}-quantile of \eqn{C_W^+}
#'       under \eqn{H_p}.}
#'     \item{\code{probs}}{The quantile level used.}
#'     \item{\code{distribution}}{Numeric vector of simulated
#'       \eqn{C_W^+} values, length \code{numsims}.}
#'     \item{\code{total_log10_lr}}{Numeric vector of simulated
#'       \eqn{\log_{10}\mathrm{LR}} totals under \eqn{H_p}.}
#'     \item{\code{numsims}}{The number of simulations actually used
#'       (profiles with non-positive total support are discarded).}
#'   }
#'
#' @details
#' The calibration uses only the \eqn{H_p} branch of
#' \code{\link{sim_lr_genetic}}: for each simulated matching profile,
#' the per-marker LR vector is extracted and passed to
#' \code{\link{concentration_index_positive}}. Simulations with no
#' positive per-marker support (empty numerator) are dropped before
#' taking the quantile.
#'
#' The default \code{probs = 0.90} matches the cutoff convention used
#' in Marsico & Egeland (in preparation): a 10\% false-flag rate under
#' \eqn{H_p} is deemed acceptable in exchange for catching concentrated
#' cases where the combined LR depends heavily on a single marker.
#'
#' @examples
#' \dontrun{
#'   library(pedtools); library(forrel)
#'   x <- linearPed(2)
#'   x <- setMarkers(x, locusAttributes = NorwegianFrequencies[1:15])
#'   x <- profileSim(x, N = 1, ids = 2)
#'   cal <- calibrate_concentration_cutoff(x, missing = 5,
#'                                         numsims = 500, probs = 0.90)
#'   cal$cutoff
#' }
#'
#' @seealso \code{\link{sim_lr_genetic}},
#'   \code{\link{concentration_index_positive}},
#'   \code{\link{fragility_report}}.
#'
#' @references
#' Marsico, F. L. & Egeland, T. (in preparation). Belief dynamics
#' during the investigative process.
#'
#' @export
calibrate_concentration_cutoff <- function(reference, missing,
                                           numsims = 1500,
                                           probs = 0.90,
                                           seed = 123,
                                           numCores = 1) {
  if (!is.numeric(probs) || length(probs) != 1 || probs <= 0 || probs >= 1)
    stop("probs must be a single numeric in (0, 1).")

  sims <- sim_lr_genetic(reference, missing = missing,
                         numsims = numsims, seed = seed,
                         numCores = numCores)

  related <- sims$Related
  cw_vals <- numeric(length(related))
  tot_vals <- numeric(length(related))
  keep <- logical(length(related))

  for (i in seq_along(related)) {
    lr_k <- related[[i]]$LRperMarker
    if (is.null(lr_k) || !is.numeric(lr_k) || length(lr_k) == 0) {
      keep[i] <- FALSE
      next
    }
    log_lrs <- log10(lr_k)
    pos <- log_lrs[log_lrs > 0]
    if (length(pos) == 0) {
      keep[i] <- FALSE
      next
    }
    cw_vals[i] <- concentration_index_positive(log_lrs)
    tot_vals[i] <- sum(log_lrs)
    keep[i] <- TRUE
  }

  cw_vals <- cw_vals[keep]
  tot_vals <- tot_vals[keep]

  if (length(cw_vals) == 0)
    stop("No simulation produced positive support under H_p; ",
         "check that the pedigree has informative markers attached.")

  list(
    cutoff          = unname(stats::quantile(cw_vals, probs = probs, na.rm = TRUE)),
    probs           = probs,
    distribution    = cw_vals,
    total_log10_lr  = tot_vals,
    numsims         = length(cw_vals)
  )
}


#' Per-case fragility report
#'
#' Produces the case-level fragility summary described in §5 of Marsico
#' & Egeland (in preparation): total weight of evidence, inclusion
#' concentration index \eqn{C_W^+}, worst-case leave-one-out loss,
#' residual support after removal of the most impactful marker, and an
#' optional comparison against a pedigree-specific cutoff (typically
#' produced by \code{\link{calibrate_concentration_cutoff}}).
#'
#' @param per_marker_lrs Named numeric vector of per-marker likelihood
#'   ratios (strictly positive). Names become the marker labels in the
#'   report.
#' @param cutoff Optional numeric scalar. If supplied, the report flags
#'   cases whose observed \eqn{C_W^+} exceeds \code{cutoff} as requiring
#'   a leave-one-out review before final reporting.
#' @param probs Optional numeric in \eqn{(0, 1)}. The quantile level at
#'   which \code{cutoff} was calibrated. Used only to enrich the
#'   reportable sentence (\code{"below the 90th-percentile cutoff"}
#'   etc.); ignored if \code{cutoff} is \code{NULL}.
#'
#' @return A list with elements
#'   \describe{
#'     \item{\code{total_log10_lr}}{Total weight of evidence
#'       \eqn{W = \sum_k \log_{10}\mathrm{LR}_k} (in bans).}
#'     \item{\code{concentration}}{Inclusion concentration index
#'       \eqn{C_W^+ \in [0, 1]}.}
#'     \item{\code{top_marker}}{Name of the single most impactful
#'       supporting marker.}
#'     \item{\code{top_log10_lr}}{Log-LR of the top marker (in bans).}
#'     \item{\code{residual_log10_lr}}{Total support that survives
#'       worst-case removal of the top marker, \eqn{(1 - C_W^+) W}.}
#'     \item{\code{flag}}{Logical. \code{TRUE} if \code{cutoff} was
#'       supplied and \eqn{C_W^+ > \text{cutoff}}; \code{NA} otherwise.}
#'     \item{\code{statement}}{Character scalar. A ready-to-paste
#'       natural-language fragility statement for the case report.}
#'   }
#'
#' @details
#' The report is a direct operationalization of Paper 1 §5:
#'
#' \itemize{
#'   \item \eqn{W} is the combined weight of evidence in bans.
#'   \item \eqn{C_W^+} answers the question ``what fraction of the
#'     inclusion support comes from a single marker?''
#'   \item \eqn{(1 - C_W^+) W} is the residual support after the
#'     single-marker worst case, i.e.\ what the analyst would be left
#'     with if a successful challenge to the top marker were sustained.
#'   \item The flag is raised when \eqn{C_W^+} exceeds the
#'     pedigree-specific cutoff, which in turn is typically the 90th
#'     percentile of the \eqn{H_p} simulation distribution produced by
#'     \code{\link{calibrate_concentration_cutoff}}.
#' }
#'
#' @examples
#' # Balanced case: 15 markers, similar contributions
#' set.seed(1)
#' lrs <- setNames(runif(15, 1.5, 3.5), paste0("M", 1:15))
#' fragility_report(lrs)
#'
#' # Concentrated case: one dominant marker
#' lrs2 <- setNames(c(1e5, rep(1.05, 14)), paste0("M", 1:15))
#' fragility_report(lrs2, cutoff = 0.20, probs = 0.90)
#'
#' @seealso \code{\link{concentration_index_positive}},
#'   \code{\link{leave_one_out}},
#'   \code{\link{calibrate_concentration_cutoff}}.
#'
#' @references
#' Marsico, F. L. & Egeland, T. (in preparation). Belief dynamics
#' during the investigative process.
#'
#' @export
fragility_report <- function(per_marker_lrs, cutoff = NULL, probs = NULL) {
  if (!is.numeric(per_marker_lrs))
    stop("per_marker_lrs must be a numeric vector.")
  if (any(per_marker_lrs <= 0)) {
    bad <- per_marker_lrs[per_marker_lrs <= 0]
    stop("All per_marker_lrs must be > 0. Values < 1 (excluding markers) are allowed; ",
         "zero or negative LR values are not. Got: ",
         paste(bad[seq_len(min(3, length(bad)))], collapse = ", "))
  }

  marker_names <- names(per_marker_lrs)
  if (is.null(marker_names))
    marker_names <- paste0("M", seq_along(per_marker_lrs))

  log_lrs <- log10(per_marker_lrs)
  W <- sum(log_lrs)
  cw <- concentration_index_positive(log_lrs)

  pos <- log_lrs > 0
  if (any(pos)) {
    top_idx <- which.max(log_lrs)
    top_name <- marker_names[top_idx]
    top_lr <- log_lrs[top_idx]
  } else {
    top_name <- NA_character_
    top_lr <- NA_real_
  }

  residual <- (1 - cw) * W

  flag <- NA
  if (!is.null(cutoff)) {
    if (!is.numeric(cutoff) || length(cutoff) != 1)
      stop("cutoff must be a single numeric value.")
    flag <- cw > cutoff
  }

  pct_label <- if (!is.null(cutoff) && !is.null(probs))
    sprintf("%.0fth-percentile", 100 * probs)
  else "calibrated"

  statement <- if (is.null(cutoff)) {
    sprintf(
      "Combined weight of evidence W = %.2f bans; C_W+ = %.3f (top marker: %s, %.2f bans); a successful challenge of the top marker would leave %.2f bans of residual support.",
      W, cw, top_name, top_lr, residual
    )
  } else if (isTRUE(flag)) {
    sprintf(
      "Combined weight of evidence W = %.2f bans; C_W+ = %.3f EXCEEDS the %s cutoff c* = %.3f (top marker: %s, %.2f bans). Leave-one-out review required: the residual support after worst-case removal is %.2f bans.",
      W, cw, pct_label, cutoff, top_name, top_lr, residual
    )
  } else {
    sprintf(
      "Combined weight of evidence W = %.2f bans; C_W+ = %.3f, below the %s cutoff c* = %.3f (top marker: %s, %.2f bans). A successful challenge of the top marker would leave %.2f bans of residual support.",
      W, cw, pct_label, cutoff, top_name, top_lr, residual
    )
  }

  list(
    total_log10_lr    = W,
    concentration     = cw,
    top_marker        = top_name,
    top_log10_lr      = unname(top_lr),
    residual_log10_lr = residual,
    flag              = flag,
    statement         = statement
  )
}
