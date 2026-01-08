#' Plot Decision Curve (FPR vs FNR)
#'
#' @description
#' Creates a scatter plot showing the trade-off between false positive rate
#' (FPR) and false negative rate (FNR) across different LR threshold values.
#' This visualization helps identify optimal decision thresholds based on
#' the relative costs of different types of errors.
#'
#' @param datasim A data.frame with columns \code{Related} and \code{Unrelated}
#'   containing LR values. Can be output from \code{\link{sim_lr_genetic}},
#'   \code{\link{sim_lr_prelim}}, \code{\link{lr_to_dataframe}}, or
#'   \code{\link{lr_combine}}.
#' @param LRmax Numeric. Maximum LR value to use as threshold. Points are
#'   generated for thresholds from 1 to LRmax. Default: 1000.
#'
#' @return A \code{ggplot2} scatter plot where:
#'   \itemize{
#'     \item X-axis: False Negative Rate (FNR) - proportion of true matches missed
#'     \item Y-axis: False Positive Rate (FPR) - proportion of non-matches incorrectly identified
#'     \item Each point represents a different LR threshold
#'   }
#'   The first and last threshold values are labeled on the plot.
#'
#' @details
#' If the input is a list (output from \code{\link{sim_lr_genetic}}), it is
#' automatically converted to a data.frame using \code{\link{lr_to_dataframe}}.
#'
#' \strong{Error Rate Definitions:}
#' \itemize{
#'   \item \emph{FPR}: Proportion of unrelated (H2) cases with LR > threshold
#'   \item \emph{FNR}: Proportion of related (H1) cases with LR < threshold
#' }
#'
#' \strong{Ideal point}: The origin (0,0) represents perfect discrimination.
#' Points closer to the origin indicate better thresholds.
#'
#' \strong{Trade-off}: Moving along the curve, decreasing FNR typically
#' increases FPR and vice versa. The optimal point depends on the relative
#' costs of false positives vs false negatives.
#'
#' @seealso
#' \code{\link{plot_lr_distribution}} for LR distribution visualization,
#' \code{\link{decision_threshold}} for computing optimal threshold,
#' \code{\link{threshold_rates}} for error rates at a specific threshold.
#'
#' @references
#' Marsico FL, Vigeland MD, Egeland T, Herrera Pinero F (2021). "Making
#' decisions in missing person identification cases with low statistical
#' power." \emph{Forensic Science International: Genetics}, 52, 102519.
#' \doi{10.1016/j.fsigen.2021.102519}
#'
#' @export
#' @import ggplot2
#'
#' @examples
#' # Using preliminary data
#' lr_sims <- sim_lr_prelim("sex", numsims = 500, seed = 123)
#' plot_decision_curve(lr_sims)
#'
#' # With lower maximum threshold for finer resolution
#' plot_decision_curve(lr_sims, LRmax = 100)

plot_decision_curve <- function(datasim, LRmax = 1000) {

  # Convert list to dataframe if needed
  if (!is.data.frame(datasim)) {
    datasim <- lr_to_dataframe(datasim)
  }

  # Avoid R CMD check notes
  x <- y <- z <- NULL

  # Extract LR values
  TPED <- datasim$Related
  RPED <- datasim$Unrelated
  nsimul <- nrow(datasim)

  # Calculate FPR and FNR for each threshold
  ValoresLR <- seq(1, LRmax, length.out = LRmax)
  FPs <- numeric(length(ValoresLR))
  FNs <- numeric(length(ValoresLR))

  for (i in seq_along(ValoresLR)) {
    FPs[i] <- sum(RPED > ValoresLR[i])
    FNs[i] <- sum(TPED < ValoresLR[i])
  }

  # Create data frame for plotting
  Datos <- data.frame(
    x = ValoresLR,
    y = FPs / nsimul,  # FPR
    z = FNs / nsimul,  # FNR
    w = FPs
  )

  # Identify first and last points for labeling
  primer_punto <- 1
  ultimo_punto <- nrow(Datos)

  # Create plot
  p <- ggplot2::ggplot(Datos, ggplot2::aes(x = z, y = y)) +
    ggplot2::geom_point(size = 2, color = "black") +
    ggplot2::geom_text(
      data = Datos[c(primer_punto, ultimo_punto), ],
      ggplot2::aes(label = x),
      hjust = 1.1, vjust = -0.5, size = 3
    ) +
    ggplot2::scale_x_continuous(name = "False Negative Rate (FNR)") +
    ggplot2::scale_y_continuous(name = "False Positive Rate (FPR)") +
    ggplot2::labs(
      title = "Decision Curve",
      subtitle = "Trade-off between FPR and FNR at different LR thresholds"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title.x = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 12),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_text(size = 10)
    )

  return(p)
}
