
#' Likelihood ratio distribution: a function for plotting expected log10(LR) distributions under relatedness and unrelatedness.
#'
#' @param datasim Input dataframe containing expected LRs for related and unrelated POIs. It should be the output from makeLRsims function.
#'
#' @export
#' @return A plot showing likelihood ratio distributions under relatedness and unrelatedness hypothesis.
#' @examples
#' library(forrel)
#' x = linearPed(2)
#' x = setMarkers(x, locusAttributes = NorwegianFrequencies[1:5])
#' x = profileSim(x, N = 1, ids = 2)
#' datasim = simLRgen(x, missing = 5, 10, 123)
#' LRdist(datasim)
#' @import dplyr
#' @import tidyr
#' @import ggplot2


LRdist = function(datasim) {
if (!is.data.frame(datasim)) {
   datasim <- simLR2dataframe(datasim)
 }

x <- y <- z <- NULL

TPED = log10(datasim$Related)
RPED = log10(datasim$Unrelated)

# Calculamos las densidades
density_TPED <- stats::density(TPED)
density_RPED <- stats::density(RPED)

# Convertimos las densidades en data frames para ggplot2
df_TPED <- data.frame(x = density_TPED$x, y = density_TPED$y)
df_RPED <- data.frame(x = density_RPED$x, y = density_RPED$y)

# Creamos el plot
ggplot2::ggplot() +
  ggplot2::geom_area(data = df_TPED, ggplot2::aes(x = x, y = y), fill = "steelblue", alpha = 0.6, color = NA) +
  ggplot2::geom_area(data = df_RPED, ggplot2::aes(x = x, y = y), fill = "#B71C1C", alpha = 0.6, color = NA) +
  ggplot2::labs(title = "LR distributions",
                x = "Log10(LR)",
                y = "Density") +
  ggplot2::theme_minimal()

}
