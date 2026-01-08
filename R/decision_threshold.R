#' Compute Optimal Decision Threshold
#'
#' @description
#' Calculates the optimal likelihood ratio (LR) threshold for classifying
#' matches versus non-matches, based on the trade-off between false positive
#' and false negative rates.
#'
#' The optimal threshold minimizes a weighted Euclidean distance that balances
#' the costs of different types of errors.
#'
#' @param datasim A data.frame with columns \code{Related} and \code{Unrelated}
#'   containing LR values. Can be output from \code{\link{sim_lr_genetic}},
#'   \code{\link{sim_lr_prelim}}, \code{\link{lr_to_dataframe}}, or
#'   \code{\link{lr_combine}}.
#' @param weight Numeric. The relative weight of false positives compared to
#'   false negatives. A value > 1 penalizes false positives more heavily.
#'   Suggested value: 10 (false positives are 10x worse than false negatives).
#'
#' @return Prints and invisibly returns the suggested LR threshold value.
#'
#' @details
#' If the input is a list (output from \code{\link{sim_lr_genetic}}), it is
#' automatically converted to a data.frame using \code{\link{lr_to_dataframe}}.
#'
#' \strong{Algorithm:}
#' The function computes the weighted Euclidean distance for each potential
#' threshold value:
#' \deqn{D = \sqrt{FNR^2 + (weight \times FPR)^2}}
#'
#' The threshold that minimizes this distance is returned as optimal.
#'
#' \strong{Weight interpretation:}
#' \itemize{
#'   \item weight = 1: Equal importance to FPR and FNR
#'   \item weight = 10: FPR is 10x more costly than FNR
#'   \item weight > 10: Very conservative (minimizes false positives)
#'   \item weight < 1: Aggressive (minimizes false negatives)
#' }
#'
#' In missing person cases, false positives (wrongly identifying someone as
#' the missing person) are typically considered more serious than false
#' negatives (failing to identify a true match), justifying weight > 1.
#'
#' @seealso
#' \code{\link{threshold_rates}} for computing error rates at a given threshold,
#' \code{\link{plot_decision_curve}} for visualizing the FPR/FNR trade-off,
#' \code{\link{plot_lr_distribution}} for LR distribution visualization.
#'
#' @references
#' Marsico FL, Vigeland MD, Egeland T, Herrera Pinero F (2021). "Making
#' decisions in missing person identification cases with low statistical
#' power." \emph{Forensic Science International: Genetics}, 52, 102519.
#' \doi{10.1016/j.fsigen.2021.102519}
#'
#' @export
#' @examples
#' # Simulate LRs
#' lr_sims <- sim_lr_prelim("sex", numsims = 500, seed = 123)
#'
#' # Find optimal threshold (FP 10x worse than FN)
#' threshold <- decision_threshold(lr_sims, weight = 10)
#'
#' # Check error rates at this threshold
#' threshold_rates(lr_sims, threshold)
#'
#' # More conservative threshold (FP 20x worse)
#' decision_threshold(lr_sims, weight = 20)

decision_threshold <- function(datasim, weight) {

  # Convert list to dataframe if needed
  if (!is.data.frame(datasim)) {
    datasim <- lr_to_dataframe(datasim)
  }

  as.data.frame(datasim)
  nsims <- nrow(datasim)
  TPED <- datasim$Related
  RPED <- datasim$Unrelated

  # Calculate FPR and FNR for each potential threshold
  ValoresLR <- seq(1, nsims, length.out = nsims)
  FPs <- numeric(10000)
  FNs <- numeric(10000)

  for (i in 1:10000) {
    FPs[i] <- sum(RPED > ValoresLR[i])
    FNs[i] <- sum(TPED < ValoresLR[i])
  }

  # Calculate weighted Euclidean distance
  Dis <- sqrt(((FNs / nsims))^2 + ((weight * FPs / nsims)^2))

  # Find threshold that minimizes distance
  Tabla <- base::data.frame(x = ValoresLR, y = Dis)
  DT <- which.min(Tabla$y)

  message(paste("Decision threshold is:", DT))
  invisible(DT)
}
