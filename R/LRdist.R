
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
#' x = profileSim(x, N = 1, ids = 2)[[1]]
#' datasim = makeLRsims(x, missing = 5, 10, 123)
#' LRdist(datasim)
#' @import plotly
#' @import dplyr
#' @import highcharter



LRdist = function(datasim) {

TPED = log10(datasim$Related)
RPED = log10(datasim$Unrelated)

hc <- hchart(
  stats::density(TPED), type = "area", 
  color = "steelblue", name = "Related"
  ) %>%
  hc_add_series(
    stats::density(RPED), type = "area",
    color = "#B71C1C", 
    name = "Unrelated"
    ) %>%
    hc_title(text = "LR distributions") %>%
    hc_xAxis(title = list(text = "Log10(LR)")) %>%
    hc_yAxis(title = list(text = "Density"))
 hc
}
