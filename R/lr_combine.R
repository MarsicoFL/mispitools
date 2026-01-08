#' Combine Likelihood Ratios from Multiple Sources
#'
#' @description
#' Combines (multiplies) likelihood ratios from two independent evidence
#' sources using the Bayesian multiplication principle. This is used to
#' integrate genetic and non-genetic evidence, or multiple non-genetic
#' variables.
#'
#' @param LRdatasim1 A data.frame with columns \code{Unrelated} and
#'   \code{Related} containing LR values from the first evidence source.
#'   Must be a data.frame (use \code{\link{lr_to_dataframe}} to convert
#'   genetic simulation output first).
#' @param LRdatasim2 A data.frame with columns \code{Unrelated} and
#'   \code{Related} containing LR values from the second evidence source.
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item \code{Unrelated}: Product of LR values under H2
#'     \item \code{Related}: Product of LR values under H1
#'   }
#'   The number of rows equals the minimum of the input data frames.
#'
#' @details
#' Under the assumption of conditional independence of evidence given
#' each hypothesis, the combined LR is the product of individual LRs:
#' \deqn{LR_{combined} = LR_1 \times LR_2}
#'
#' This follows from Bayes' theorem and is valid when the evidence sources
#' are conditionally independent given the hypothesis.
#'
#' \strong{Important:} Both inputs must be data.frames with the same
#' structure. If using output from \code{\link{sim_lr_genetic}}, first
#' convert it using \code{\link{lr_to_dataframe}}.
#'
#' @seealso
#' \code{\link{sim_lr_genetic}} for genetic LR simulations,
#' \code{\link{sim_lr_prelim}} for non-genetic LR simulations,
#' \code{\link{lr_to_dataframe}} for converting genetic simulations,
#' \code{\link{plot_lr_distribution}} for visualizing combined distributions.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @examples
#' # Simulate LRs from two different variables
#' lr_sex <- sim_lr_prelim("sex", numsims = 500, seed = 123)
#' lr_age <- sim_lr_prelim("age", numsims = 500, seed = 456)
#'
#' # Combine the evidence
#' lr_combined <- lr_combine(lr_sex, lr_age)
#' head(lr_combined)
#'
#' # Compare distributions
#' summary(log10(lr_sex$Related))
#' summary(log10(lr_combined$Related))
#'
#' # Visualize combined distribution
#' plot_lr_distribution(lr_combined)
#'
#' # Combining genetic and non-genetic evidence
#' library(forrel)
#' x <- linearPed(2)
#' x <- setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x <- profileSim(x, N = 1, ids = 2)
#'
#' # Simulate genetic LRs and convert to dataframe
#' lr_genetic <- sim_lr_genetic(x, missing = 5, numsims = 100, seed = 123)
#' lr_genetic_df <- lr_to_dataframe(lr_genetic)
#'
#' # Simulate non-genetic LRs
#' lr_prelim <- sim_lr_prelim("sex", numsims = 100, seed = 123)
#'
#' # Combine both sources
#' lr_total <- lr_combine(lr_genetic_df, lr_prelim)

lr_combine <- function(LRdatasim1, LRdatasim2) {

  # Ensure both inputs are data.frames
  if (!is.data.frame(LRdatasim1)) {
    stop("LRdatasim1 must be a data.frame. Use lr_to_dataframe() to convert genetic simulations.")
  }
  if (!is.data.frame(LRdatasim2)) {
    stop("LRdatasim2 must be a data.frame. Use lr_to_dataframe() to convert genetic simulations.")
  }

  # Check required columns
  required_cols <- c("Unrelated", "Related")
  if (!all(required_cols %in% names(LRdatasim1))) {
    stop("LRdatasim1 must have columns 'Unrelated' and 'Related'")
  }
  if (!all(required_cols %in% names(LRdatasim2))) {
    stop("LRdatasim2 must have columns 'Unrelated' and 'Related'")
  }

  # Check row counts match
  n1 <- nrow(LRdatasim1)
  n2 <- nrow(LRdatasim2)
  if (n1 != n2) {
    stop("LRdatasim1 and LRdatasim2 must have the same number of rows (",
         n1, " vs ", n2, "). ",
         "LR combination requires paired simulations under the same scenarios.")
  }

  # Multiply LRs (element-wise)
  # This assumes conditional independence of evidence given each hypothesis
  a <- as.numeric(unlist(LRdatasim1$Unrelated)) * as.numeric(unlist(LRdatasim2$Unrelated))
  b <- as.numeric(unlist(LRdatasim1$Related)) * as.numeric(unlist(LRdatasim2$Related))

  # Combine into dataframe
  LRsimulated <- cbind(a, b)
  colnames(LRsimulated) <- c("Unrelated", "Related")

  structure(as.data.frame(LRsimulated))
}
