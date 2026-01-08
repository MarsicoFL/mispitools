#' Simulate LR Distributions for Pigmentation Traits
#'
#' @description
#' Simulates likelihood ratio (LR) distributions for combined pigmentation
#' traits (hair, skin, and eye color) under both hypotheses. Uses pre-computed
#' LRs from \code{\link{lr_compute_pigmentation}}.
#'
#' @param df A data.frame with columns \code{numerators}, \code{f_h_s_y},
#'   and \code{LR}. Typically output from \code{\link{lr_compute_pigmentation}}.
#' @param seed Integer. Random seed for reproducibility. Default: 1234.
#' @param nsim Integer. Number of LR values to simulate per hypothesis.
#'   Default: 500.
#'
#' @return A data.frame with two columns:
#'   \itemize{
#'     \item \code{Unrelated}: LR values simulated under H2 (sampling proportional
#'       to population frequencies)
#'     \item \code{Related}: LR values simulated under H1 (sampling proportional
#'       to conditioned probabilities)
#'   }
#'
#' @details
#' The function samples LR values with probabilities proportional to:
#' \itemize{
#'   \item \emph{H2 (Unrelated)}: Population frequencies (\code{f_h_s_y})
#'   \item \emph{H1 (Related)}: Conditioned probabilities (\code{numerators})
#' }
#'
#' This simulates the expected distribution of LRs when comparing the MP's
#' traits against either random individuals (H2) or the true match (H1).
#'
#' @seealso
#' \code{\link{sim_reference_pop}} for generating population data,
#' \code{\link{lr_compute_pigmentation}} for computing input LRs,
#' \code{\link{plot_lr_distribution}} for visualization.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @examples
#' # Full workflow for pigmentation LRs
#' pop_data <- sim_reference_pop(n = 500, seed = 123)
#' conditioned <- compute_conditioned_prop(pop_data, 1, 1, 1, 0.01, 0.01, 0.01)
#' unconditioned <- compute_reference_prop(pop_data)
#' lrs <- lr_compute_pigmentation(conditioned, unconditioned)
#'
#' # Simulate LR distribution
#' lr_dist <- lr_pigmentation(lrs, nsim = 500, seed = 456)
#' head(lr_dist)
#'
#' # Visualize
#' plot_lr_distribution(lr_dist)

lr_pigmentation <- function(df, seed = 1234, nsim = 500) {

  set.seed(seed)

  LR <- df$LR
  f_h_s_y <- df$f_h_s_y
  numerators <- df$numerators

  # Normalize probabilities
  prob_unrelated <- f_h_s_y / sum(f_h_s_y)
  prob_related <- numerators / sum(numerators)

  # Sample LRs proportional to probabilities
  Unrelated <- sample(LR, nsim, replace = TRUE, prob = prob_unrelated)
  Related <- sample(LR, nsim, replace = TRUE, prob = prob_related)

  result_df <- data.frame(Unrelated = Unrelated, Related = Related)

  return(result_df)
}
