#' Compute Error Rates at a Specific Threshold
#'
#' @description
#' Calculates error rates and performance metrics for a given likelihood ratio
#' (LR) threshold, including:
#' \itemize{
#'   \item False Positive Rate (FPR)
#'   \item False Negative Rate (FNR)
#'   \item Matthews Correlation Coefficient (MCC)
#' }
#'
#' @param datasim A data.frame with columns \code{Related} and \code{Unrelated}
#'   containing LR values. Can be output from \code{\link{sim_lr_genetic}},
#'   \code{\link{sim_lr_prelim}}, \code{\link{lr_to_dataframe}}, or
#'   \code{\link{lr_combine}}.
#' @param threshold Numeric. The LR threshold value for which to compute
#'   error rates. Cases with LR > threshold are classified as matches.
#'
#' @return Prints the error rates and MCC, and invisibly returns a named list
#'   with components:
#'   \itemize{
#'     \item \code{FNR}: False Negative Rate
#'     \item \code{FPR}: False Positive Rate
#'     \item \code{TPR}: True Positive Rate
#'     \item \code{TNR}: True Negative Rate
#'     \item \code{MCC}: Matthews Correlation Coefficient
#'   }
#'
#' @details
#' If the input is a list (output from \code{\link{sim_lr_genetic}}), it is
#' automatically converted to a data.frame using \code{\link{lr_to_dataframe}}.
#'
#' \strong{Metrics:}
#' \itemize{
#'   \item \emph{FPR}: Proportion of unrelated cases incorrectly classified as matches
#'     (LR > threshold when H2 is true)
#'   \item \emph{FNR}: Proportion of related cases incorrectly classified as non-matches
#'     (LR < threshold when H1 is true)
#'   \item \emph{TPR}: 1 - FNR (sensitivity, recall)
#'   \item \emph{TNR}: 1 - FPR (specificity)
#'   \item \emph{MCC}: Matthews Correlation Coefficient, ranges from -1 to +1:
#'     \itemize{
#'       \item +1: Perfect classification
#'       \item 0: Random classification
#'       \item -1: Completely wrong classification
#'     }
#' }
#'
#' @seealso
#' \code{\link{decision_threshold}} for finding optimal threshold,
#' \code{\link{plot_decision_curve}} for visualizing the FPR/FNR trade-off.
#'
#' @references
#' Marsico FL, Vigeland MD, Egeland T, Herrera Pinero F (2021). "Making
#' decisions in missing person identification cases with low statistical
#' power." \emph{Forensic Science International: Genetics}, 52, 102519.
#' \doi{10.1016/j.fsigen.2021.102519}
#'
#' Matthews BW (1975). "Comparison of the predicted and observed secondary
#' structure of T4 phage lysozyme." \emph{Biochimica et Biophysica Acta},
#' 405(2), 442-451.
#'
#' @export
#' @examples
#' # Simulate LRs
#' lr_sims <- sim_lr_prelim("sex", numsims = 500, seed = 123)
#'
#' # Check error rates at threshold = 10
#' rates <- threshold_rates(lr_sims, threshold = 10)
#'
#' # Access individual metrics
#' rates$FPR
#' rates$MCC
#'
#' # Compare different thresholds
#' threshold_rates(lr_sims, threshold = 5)
#' threshold_rates(lr_sims, threshold = 50)
#' threshold_rates(lr_sims, threshold = 100)

threshold_rates <- function(datasim, threshold) {

  # Input validation
  if (missing(threshold) || !is.numeric(threshold) || length(threshold) != 1) {
    stop("threshold must be a single numeric value")
  }

  # Convert list to dataframe if needed
  if (!is.data.frame(datasim)) {
    datasim <- lr_to_dataframe(datasim)
  }

  # Validate required columns
  if (!all(c("Related", "Unrelated") %in% names(datasim))) {
    stop("datasim must have columns 'Related' and 'Unrelated'")
  }

  nsims <- nrow(datasim)

  if (nsims < 1) {
    stop("datasim must have at least 1 row")
  }

  TPED <- datasim$Related
  RPED <- datasim$Unrelated

  # Calculate counts (as numeric to prevent integer overflow in MCC)
  TP <- as.numeric(sum(TPED > threshold))
  TN <- as.numeric(sum(RPED < threshold))
  FP <- as.numeric(sum(RPED > threshold))
  FN <- as.numeric(sum(TPED < threshold))

  # Calculate rates
  FPR <- FP / nsims
  FNR <- FN / nsims
  TPR <- TP / nsims
  TNR <- TN / nsims

  # Matthews Correlation Coefficient using counts
  # MCC = (TP*TN - FP*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))
  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  MCC <- if (denom == 0) 0 else (TP * TN - FP * FN) / denom

  # Print results
  message(paste("FNR =", round(FNR, 4),
                ";  FPR =", round(FPR, 4),
                ";  MCC =", round(MCC, 4)))

  # Return invisibly
  invisible(list(
    FNR = FNR,
    FPR = FPR,
    TPR = TPR,
    TNR = TNR,
    MCC = MCC
  ))
}
