#' Plot Likelihood Ratio Distributions
#'
#' @description
#' Creates a density plot showing the expected log10(LR) distributions under
#' both hypotheses:
#' \itemize{
#'   \item H1 (Related/Blue): Distribution when POI is the missing person
#'   \item H2 (Unrelated/Red): Distribution when POI is unrelated
#' }
#'
#' This visualization helps assess the discriminatory power of the evidence
#' and identify potential overlap between the two hypotheses.
#'
#' @param datasim A data.frame with columns \code{Related} and \code{Unrelated}
#'   containing LR values. Can be output from \code{\link{sim_lr_genetic}},
#'   \code{\link{sim_lr_prelim}}, \code{\link{lr_to_dataframe}}, or
#'   \code{\link{lr_combine}}.
#'
#' @return A \code{ggplot2} object showing overlaid density curves.
#'   Blue area represents H1 (Related), red area represents H2 (Unrelated).
#'
#' @details
#' If the input is a list (output from \code{\link{sim_lr_genetic}}), it is
#' automatically converted to a data.frame using \code{\link{lr_to_dataframe}}.
#'
#' The x-axis shows log10(LR), which is more interpretable than raw LR values:
#' \itemize{
#'   \item log10(LR) = 0 means LR = 1 (neutral evidence)
#'   \item log10(LR) > 0 means evidence favors H1 (related)
#'   \item log10(LR) < 0 means evidence favors H2 (unrelated)
#' }
#'
#' Less overlap between distributions indicates better discrimination.
#'
#' @seealso
#' \code{\link{plot_decision_curve}} for FPR/FNR trade-off visualization,
#' \code{\link{decision_threshold}} for computing optimal thresholds,
#' \code{\link{sim_lr_genetic}}, \code{\link{sim_lr_prelim}} for generating input.
#'
#' @references
#' Marsico FL, Vigeland MD, Egeland T, Herrera Pinero F (2021). "Making
#' decisions in missing person identification cases with low statistical
#' power." \emph{Forensic Science International: Genetics}, 52, 102519.
#' \doi{10.1016/j.fsigen.2021.102519}
#'
#' @export
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' # Using preliminary data
#' lr_sims <- sim_lr_prelim("sex", numsims = 500, seed = 123)
#' plot_lr_distribution(lr_sims)
#'
#' # Using genetic data
#' library(forrel)
#' x <- linearPed(2)
#' x <- setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x <- profileSim(x, N = 1, ids = 2)
#' lr_genetic <- sim_lr_genetic(x, missing = 5, numsims = 50, seed = 123)
#' plot_lr_distribution(lr_genetic)

plot_lr_distribution <- function(datasim) {

  # Convert list to dataframe if needed
  if (!is.data.frame(datasim)) {
    datasim <- lr_to_dataframe(datasim)
  }

  # Avoid R CMD check notes
  x <- y <- NULL

  # Transform to log10 scale
  TPED <- log10(datasim$Related)
  RPED <- log10(datasim$Unrelated)

  # Calculate kernel density estimates
  density_TPED <- stats::density(TPED)
  density_RPED <- stats::density(RPED)

  # Convert to data frames for ggplot2
  df_TPED <- data.frame(x = density_TPED$x, y = density_TPED$y)
  df_RPED <- data.frame(x = density_RPED$x, y = density_RPED$y)

  # Create plot
  ggplot2::ggplot() +
    ggplot2::geom_area(data = df_TPED,
                       ggplot2::aes(x = x, y = y),
                       fill = "steelblue", alpha = 0.6, color = NA) +
    ggplot2::geom_area(data = df_RPED,
                       ggplot2::aes(x = x, y = y),
                       fill = "#B71C1C", alpha = 0.6, color = NA) +
    ggplot2::labs(
      title = "Likelihood Ratio Distributions",
      subtitle = "Blue = H1 (Related) | Red = H2 (Unrelated)",
      x = expression(log[10](LR)),
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(size = 11)
    )
}
