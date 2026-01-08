#' Plot Conditional Probability Tables Comparison
#'
#' @description
#' Creates a three-panel visualization comparing conditional probability
#' tables (CPTs) and their resulting likelihood ratios:
#' \itemize{
#'   \item Panel A: P(D|H2) - Population-based probabilities
#'   \item Panel B: P(D|H1) - Missing person-based probabilities
#'   \item Panel C: log10(LR) - Likelihood ratios for each combination
#' }
#'
#' This visualization helps understand how different combinations of
#' sex, age group, and hair color contribute to the likelihood ratio.
#'
#' @param CPT_POP Matrix. Population-based conditional probability table,
#'   typically output from \code{\link{cpt_population}}.
#' @param CPT_MP Matrix. Missing person-based conditional probability table,
#'   typically output from \code{\link{cpt_missing_person}}.
#'
#' @return A \code{ggplot2} object with three panels arranged horizontally,
#'   showing heatmaps with cell values annotated.
#'
#' @details
#' The heatmaps use a blue gradient where darker colors indicate higher
#' values (higher probabilities or higher LRs).
#'
#' Each cell is labeled with its value rounded to 2 decimal places.
#'
#' The LR panel (C) shows log10(LR), where:
#' \itemize{
#'   \item Positive values (blue) favor H1 (related)
#'   \item Negative values favor H2 (unrelated)
#'   \item Zero indicates neutral evidence
#' }
#'
#' Row labels indicate sex and age group combinations:
#' \itemize{
#'   \item F-T1: Female, age within MP range
#'   \item F-T0: Female, age outside MP range
#'   \item M-T1: Male, age within MP range
#'   \item M-T0: Male, age outside MP range
#' }
#'
#' Column labels indicate hair color categories (1-5).
#'
#' @seealso
#' \code{\link{cpt_population}} for creating the H2 table,
#' \code{\link{cpt_missing_person}} for creating the H1 table.
#'
#' @references
#' Marsico FL, et al. (2023). "Likelihood ratios for non-genetic evidence
#' in missing person cases." \emph{Forensic Science International: Genetics},
#' 66, 102891. \doi{10.1016/j.fsigen.2023.102891}
#'
#' @export
#' @import reshape2
#' @import patchwork
#' @import ggplot2
#' @importFrom graphics par
#' @examples
#' # Create both CPTs
#' cpt_h2 <- cpt_population()
#' cpt_h1 <- cpt_missing_person(MPs = "F", MPc = 1)
#'
#' # Visualize comparison
#' plot_cpt(cpt_h2, cpt_h1)
#'
#' # Different MP characteristics
#' cpt_h1_male <- cpt_missing_person(MPs = "M", MPc = 3)
#' plot_cpt(cpt_h2, cpt_h1_male)

plot_cpt <- function(CPT_POP, CPT_MP) {

  # Avoid R CMD check notes
  Var1 <- Var2 <- value <- NULL

  graphics::par(mfrow = c(2, 1), mar = c(2, 4, 4, 2))

  # Panel A: Population probabilities (H2)
  POP <- reshape2::melt(CPT_POP)
  p1 <- ggplot2::ggplot(POP, ggplot2::aes(x = Var2, y = Var1)) +
    ggplot2::geom_raster(ggplot2::aes(fill = value)) +
    ggplot2::scale_fill_gradient(low = "grey90", high = "steelblue", limits = c(0, 1)) +
    ggplot2::labs(
      x = "Hair color",
      y = "Sex - Age group",
      title = "P(D|H2)",
      fill = "Prob"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 11, angle = 0, vjust = 0.3),
      axis.text.y = ggplot2::element_text(size = 11),
      plot.title = ggplot2::element_text(size = 13, face = "bold")
    ) +
    ggplot2::geom_text(ggplot2::aes(label = format(round(value, 2), nsmall = 2)), size = 3)

  # Panel B: Missing person probabilities (H1)
  MP <- reshape2::melt(CPT_MP)
  p2 <- ggplot2::ggplot(MP, ggplot2::aes(x = Var2, y = Var1)) +
    ggplot2::geom_raster(ggplot2::aes(fill = value)) +
    ggplot2::scale_fill_gradient(low = "grey90", high = "steelblue", limits = c(0, 1)) +
    ggplot2::labs(
      x = "Hair color",
      y = "Sex - Age group",
      title = "P(D|H1)",
      fill = "Prob"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 11, angle = 0, vjust = 0.3),
      axis.text.y = ggplot2::element_text(size = 11),
      plot.title = ggplot2::element_text(size = 13, face = "bold")
    ) +
    ggplot2::geom_text(ggplot2::aes(label = format(round(value, 2), nsmall = 2)), size = 3)

  # Panel C: Likelihood ratios
  LRtable <- log10(CPT_MP / CPT_POP)
  LR <- reshape2::melt(LRtable)
  p3 <- ggplot2::ggplot(LR, ggplot2::aes(x = Var2, y = Var1)) +
    ggplot2::geom_raster(ggplot2::aes(fill = value)) +
    ggplot2::scale_fill_gradient(low = "grey90", high = "steelblue") +
    ggplot2::labs(
      x = "Hair color",
      y = "Sex - Age group",
      title = expression(log[10](LR)),
      fill = "log10(LR)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 11, angle = 0, vjust = 0.3),
      axis.text.y = ggplot2::element_text(size = 11),
      plot.title = ggplot2::element_text(size = 13, face = "bold")
    ) +
    ggplot2::geom_text(ggplot2::aes(label = format(round(value, 2), nsmall = 2)), size = 3)

  # Combine panels
  p <- (p1 + p2 + p3) + patchwork::plot_annotation(tag_levels = 'A')

  return(p)
}
