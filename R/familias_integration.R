#' Belief trajectory metrics from a \code{Familias::FamiliasPosterior} result
#'
#' Convenience wrapper that extracts the per-marker likelihood ratios from
#' the result list returned by \code{Familias::FamiliasPosterior} and
#' computes the full belief-trajectory machinery of Marsico & Egeland
#' (in preparation): the binary belief trajectory, the trajectory metrics
#' (entropy, per-step KL divergence, total-variation path length, and the
#' three concentration measures), the signed concentration index
#' \eqn{C_W^+} for fragility of inclusions, and the per-marker
#' leave-one-out analysis.
#'
#' @param familias_result A list with the structure returned by
#'   \code{Familias::FamiliasPosterior}. Must contain a numeric matrix named
#'   \code{LRperMarker} with rows indexing markers and columns indexing
#'   candidate pedigrees. Column \code{ref_pedigree} is assumed to contain
#'   the reference hypothesis (typically all ones after normalization).
#' @param test_pedigree Integer or character. Index or name of the
#'   alternative pedigree whose per-marker LRs (as a ratio against the
#'   reference pedigree) should be used to build the belief trajectory.
#'   Defaults to 2 (i.e., the first non-reference pedigree), which is the
#'   standard arrangement in two-pedigree comparisons.
#' @param ref_pedigree Integer or character. Index or name of the
#'   reference pedigree. Defaults to 1. Used for sanity checks only; the
#'   actual LR-per-marker values for \code{test_pedigree} are taken from
#'   \code{familias_result$LRperMarker} directly because
#'   \code{FamiliasPosterior} already normalizes the matrix against the
#'   reference column.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{lrs}}{Named numeric vector of per-marker likelihood
#'     ratios for the chosen test pedigree against the reference.}
#'   \item{\code{trajectory}}{Data frame from
#'     \code{\link{binary_belief_trajectory}} with the posterior at each
#'     step, cumulative log-LR, and the per-step log-LR.}
#'   \item{\code{metrics}}{List from \code{\link{trajectory_metrics}}
#'     with entropy, per-step KL divergence, cumulative KL from prior,
#'     per-step total-variation, path length, and the three
#'     concentration measures.}
#'   \item{\code{concentration_positive}}{Numeric scalar. The signed
#'     concentration index \eqn{C_W^+} restricted to positive
#'     (supporting) per-marker contributions; the primary fragility
#'     diagnostic of the Belief Dynamics framework.}
#'   \item{\code{leave_one_out}}{Data frame from
#'     \code{\link{leave_one_out}} giving the per-marker fragility table.}
#' }
#'
#' @details
#' \code{Familias::FamiliasPosterior} computes posterior probabilities of
#' candidate pedigrees given DNA evidence. Its return value includes a
#' \code{LRperMarker} matrix of per-locus likelihood ratios already
#' normalized against the reference pedigree. This wrapper simply selects
#' the relevant column and hands the resulting per-marker LR vector to
#' the Belief Dynamics trajectory machinery, so that a user of Familias
#' can obtain trajectory metrics and fragility diagnostics in a single
#' function call without manually constructing the per-marker sequence.
#'
#' If \code{Familias} is not installed, this function does not depend on
#' it --- the input is expected to be a list with the structure described
#' above, which can equally well be constructed by hand for testing or
#' from other sources (e.g., \code{forrel::missingPersonLR}).
#'
#' @examples
#' # Synthetic Familias-like result with 6 markers and 2 pedigrees
#' fam <- list(
#'   LRperMarker = matrix(
#'     c(1, 1, 1, 1, 1, 1,          # reference pedigree, col 1
#'       5.2, 1.8, 12.0, 3.1, 0.8, 2.7), # test pedigree,     col 2
#'     nrow = 6, ncol = 2,
#'     dimnames = list(c("D3S1358", "TH01", "D21S11",
#'                       "D18S51",  "CSF", "vWA"),
#'                     c("Unrelated", "GrandparentGrandchild"))
#'   )
#' )
#' result <- familias_trajectory(fam, test_pedigree = 2)
#' result$concentration_positive
#' result$leave_one_out
#'
#' @seealso \code{\link{binary_belief_trajectory}},
#'   \code{\link{trajectory_metrics}},
#'   \code{\link{concentration_index_positive}},
#'   \code{\link{leave_one_out}}.
#'
#' @references
#' Marsico, F. L. & Egeland, T. (in preparation).
#' Belief dynamics during the investigative process.
#' Egeland, T., Mostad, P. & Simonsson, I. (2015).
#' Relationship inference with Familias and R. Academic Press.
#'
#' @export
familias_trajectory <- function(familias_result,
                                test_pedigree = 2,
                                ref_pedigree = 1) {
  if (!is.list(familias_result))
    stop("familias_result must be a list (as returned by Familias::FamiliasPosterior)")
  if (is.null(familias_result$LRperMarker))
    stop("familias_result must contain a matrix named 'LRperMarker'")
  lrm <- familias_result$LRperMarker
  if (!is.matrix(lrm) || !is.numeric(lrm))
    stop("familias_result$LRperMarker must be a numeric matrix")
  if (ncol(lrm) < 2)
    stop("familias_result$LRperMarker must have at least 2 columns (reference + at least one test pedigree)")

  # Resolve test_pedigree by name or index
  if (is.character(test_pedigree)) {
    idx <- match(test_pedigree, colnames(lrm))
    if (is.na(idx))
      stop(sprintf("test_pedigree '%s' not found in column names of LRperMarker",
                   test_pedigree))
    test_pedigree <- idx
  }
  if (test_pedigree < 1 || test_pedigree > ncol(lrm))
    stop("test_pedigree index out of range")
  if (is.character(ref_pedigree)) {
    idx <- match(ref_pedigree, colnames(lrm))
    if (is.na(idx))
      stop(sprintf("ref_pedigree '%s' not found in column names of LRperMarker",
                   ref_pedigree))
    ref_pedigree <- idx
  }
  if (test_pedigree == ref_pedigree)
    stop("test_pedigree and ref_pedigree must be different")

  lrs <- lrm[, test_pedigree]
  if (!is.null(rownames(lrm)))
    names(lrs) <- rownames(lrm)

  if (any(lrs <= 0))
    stop("Per-marker LRs must be strictly positive; check that Familias returned a proper LRperMarker matrix")

  # Binary trajectory (data frame with posterior at each step)
  traj_df <- binary_belief_trajectory(lrs)

  # General-form trajectory for the metrics machinery
  prior <- c(0.5, 0.5)
  lr_list <- lapply(unname(lrs), function(r) c(r, 1))
  traj_mat <- belief_trajectory(prior, lr_list)
  metrics <- trajectory_metrics(traj_mat)

  list(
    lrs = lrs,
    trajectory = traj_df,
    metrics = metrics,
    concentration_positive = concentration_index_positive(log10(lrs)),
    leave_one_out = leave_one_out(lrs)
  )
}
